#include <filesystem>
#include <string>

#include <m_pd.h>

#include <m_imp.h>

// libspatialaudio
#include <AmbisonicBinauralizer.h>
#include <AmbisonicDecoder.h>
#include <AmbisonicEncoderDist.h>
#include <BFormat.h>

namespace fs = std::filesystem;

static t_class *ambi_class;

// ==============================================
class ambi_tilde {
  public:
    t_object xObj;
    t_canvas *xCanvas;
    std::string patchDir;
    std::string exterDir;
    t_sample sample;
    unsigned blockSize;
    unsigned sampleRate;
    unsigned tailLength;

    // multichannels
    unsigned Chs;
    bool multichannel;

    //
    unsigned nChIn;
    unsigned nChOut;
    std::string SpeakerConfig;
    Amblib_SpeakerSetUps SpeakerSetUp;
    bool NeedSpeakersPosition = false;
    bool SpeakerPositionSet = false;

    // Position
    PolarPoint Position;

    // Encoder config
    CAmbisonicEncoder *Encoder;

    // Decoder config
    CAmbisonicDecoder *Decoder;

    // Decoder config
    CBFormat *BFormat;

    // Binaural config
    bool binauralSet = false;
    CAmbisonicBinauralizer *Binauralizer;
    std::string hrtfPath;
    bool lowcpu;

    //
    t_sample **SpeakersArray;
    t_sample **SpeakersSetup;

    bool DecoderConfigured = false;
    bool binaural = false;

    t_outlet **outlets;
};

// ─────────────────────────────────────
static void ambi_configure(ambi_tilde *x, int blockSize, int sampleRate) {
    if (blockSize == 0 || sampleRate == 0) {
        pd_error(x, "[ambi~]: Init Dsp first!");
        return;
    }
    x->DecoderConfigured = false;

    int oldSpeakersCount = x->Decoder->GetSpeakerCount();

    float fadeMs = 1000.f * blockSize / sampleRate;
    bool sucess = x->BFormat->Configure(1, true, blockSize);
    if (!sucess) {
        pd_error(x, "[ambi~]: Error configuring BFormat");
        return;
    }

    sucess = x->Encoder->Configure(1, true, sys_getsr(), fadeMs);
    if (!sucess) {
        pd_error(x, "[ambi~]: Error configuring Encoder");
        return;
    }

    unsigned int DecoderSpeakers;
    if (x->binaural) {
        sucess = x->Binauralizer->Configure(1, true, sys_getsr(), blockSize,
                                            x->tailLength, x->hrtfPath);
        if (!sucess) {
            pd_error(x, "[ambi~]: Error configuring Binauralizer");
            return;
        }
    } else {
        sucess = x->Decoder->Configure(1, true, blockSize, sys_getsr(),
                                       x->SpeakerSetUp, x->nChOut);
        if (!sucess) {
            pd_error(x, "[ambi~]: Error configuring Decoder");
            return;
        }
    }
    DecoderSpeakers = x->Decoder->GetSpeakerCount();

    if (x->SpeakersArray) {
        for (int niSpeaker = 0; niSpeaker < oldSpeakersCount; niSpeaker++) {
            delete[] x->SpeakersArray[niSpeaker];
        }
        delete[] x->SpeakersArray;
    }

    if (x->NeedSpeakersPosition && x->SpeakersSetup) {
        for (int i = 0; i < x->nChOut; i++) {
            delete[] x->SpeakersSetup[i];
        }
        x->SpeakersSetup = new t_sample *[x->nChOut];
        for (int i = 0; i < x->nChOut; i++) {
            x->SpeakersSetup[i] = new t_sample[2];
            x->SpeakersSetup[i][0] = i;
            x->SpeakersSetup[i][1] = 0;
        }
    }
    printf("position ok\n");

    x->SpeakersArray = new t_sample *[DecoderSpeakers];
    for (int niSpeaker = 0; niSpeaker < DecoderSpeakers; niSpeaker++) {
        x->SpeakersArray[niSpeaker] = new t_sample[blockSize];
    }
    x->DecoderConfigured = true;
    printf("ok\n");
}

// ==============================================
// returns true on sucess
static bool ambi_speakerconfig(ambi_tilde *x) {

    if (x->nChOut == 2) {
        x->SpeakerConfig = "stereo";
    } else if (x->nChOut == 4) {
        x->SpeakerConfig = "quad";
    } else {
        pd_error(x, "[ambi~]: You are using a custom speaker configuration "
                    "but you haven't specified the number of speakers. Use "
                    "[-s custom] to set.");
        return false;
    }

    std::string config = x->SpeakerConfig;
    if (config == "stereo") {
        x->SpeakerSetUp = Amblib_SpeakerSetUps::kAmblib_Stereo;
        x->SpeakerPositionSet = true;
        return true;
    } else if (config == "quad") {
        x->SpeakerSetUp = Amblib_SpeakerSetUps::kAmblib_Quad;
        x->SpeakerPositionSet = true;
        return true;
    } else if (config == "5.1") {
        x->SpeakerSetUp = Amblib_SpeakerSetUps::kAmblib_51;
        x->SpeakerPositionSet = true;
        return true;
    } else if (config == "7.1") {
        x->SpeakerSetUp = Amblib_SpeakerSetUps::kAmblib_71;
        x->SpeakerPositionSet = true;
        return true;
    } else if (config == "custom") {
        x->SpeakerSetUp = Amblib_SpeakerSetUps::kAmblib_CustomSpeakerSetUp;
        x->NeedSpeakersPosition = true;
        x->SpeakerPositionSet = false;
        return true;
    } else {
        pd_error(x, "[ambi~]: invalid speaker configuration");
        x->SpeakerPositionSet = false;
        return false;
    }
}

// ─────────────────────────────────────
static void ambi_binauralconfig(ambi_tilde *x, t_symbol *s, int argc,
                                t_atom *argv) {
    if (argc == 2) {
        bool binaural = atom_getfloat(argv + 1);
        if (binaural) {
            x->binaural = true;
            if (!x->binauralSet) {
                x->Binauralizer->Configure(1, true, sys_getsr(), x->blockSize,
                                           x->tailLength, x->hrtfPath);
            }
            post("[ambi~]: Binaural on");
            x->binauralSet = true;
        } else {
            x->binaural = false;
            post("[ambi~]: Binaural off");
        }
    }

    if (argc == 3) {
        std::string method = atom_getsymbol(argv + 1)->s_name;
        if (method == "hrtf") {
            if (argc < 3) {
                pd_error(x, "[ambi~]: HRTF file path not provided");
                return;
            }
            std::string file =
                x->patchDir + "/" + atom_getsymbol(argv + 2)->s_name;
            if (fs::exists(file)) {
                x->hrtfPath = file;
            } else {
                pd_error(x, "[ambi~]: HRTF file not found");
                return;
            }
        } else if (method == "lowcpu") {
            bool lowcpu = atom_getfloat(argv + 2);
            if (lowcpu) {
                x->lowcpu = true;
                post("[ambi~]: Low CPU mode on");
            } else {
                x->lowcpu = false;
                post("[ambi~]: Low CPU mode off");
            }
        }

        if (x->binaural) {
            x->Binauralizer->Configure(1, true, sys_getsr(), x->blockSize,
                                       x->tailLength, x->hrtfPath);
        }
    }
}

// ─────────────────────────────────────
static void ambi_encoderconfig(ambi_tilde *x, t_symbol *s, int argc,
                               t_atom *argv) {
    std::string method = atom_getsymbol(argv + 1)->s_name;
    if (method == "gain") {
        if (argc < 3) {
            pd_error(x, "[ambi~]: Gain value not provided");
            return;
        }
        float gain = atom_getfloat(argv + 2);
        x->Encoder->SetGain(gain);
    } else if (method == "orderweight") {
        if (argc < 4) {
            pd_error(x, "[ambi~]: Order and weight not provided");
            return;
        }
        unsigned order = atom_getint(argv + 2);
        float weight = atom_getfloat(argv + 3);
        post("[ambi~]: Order %d Weight %f", order, weight);
        x->Encoder->SetOrderWeight(order, weight);
    } else if (method == "coeff") {
        if (argc < 4) {
            pd_error(x, "[ambi~]: Channel and Coeff not provided");
            return;
        }
        unsigned channel = atom_getint(argv + 2);
        float coeff = atom_getfloat(argv + 3);
        post("[ambi~]: Channel %d Coeff %f", channel, coeff);
        x->Encoder->SetCoefficient(channel, coeff);
    } else {
        pd_error(x, "[ambi~]: invalid method, options are: gain, orderweight, "
                    "coeff");
        return;
    }

    x->Encoder->Refresh();
}

// ─────────────────────────────────────
static void ambi_decoderconfig(ambi_tilde *x, t_symbol *s, int argc,
                               t_atom *argv) {

    if (!x->DecoderConfigured) {
        pd_error(x, "[ambi~]: Decoder not configured yet, turn the DSP on");
    }

    std::string method = atom_getsymbol(argv + 1)->s_name;
    if (method == "orderweight") {
        if (argc < 5) {
            pd_error(x, "[ambi~]: Order and weight not provided");
            return;
        }
        unsigned speaker = atom_getint(argv + 2) - 1;
        unsigned order = atom_getint(argv + 3);
        float weight = atom_getfloat(argv + 4);
        if (speaker > x->Decoder->GetSpeakerCount() - 1) {
            pd_error(x,
                     "[ambi~]: There is just %d speakers in the current "
                     "configuration",
                     x->Decoder->GetSpeakerCount());
            return;
        } else if (speaker < 0) {
            pd_error(x, "[ambi~]: Speaker must be a positive number");
            return;
        }

        if (order < 1) {
            pd_error(x, "[ambi~]: Order must be from 1 to 7");
            return;
        }

        post("[ambi~]: Speaker %d Order %d Weight %f", speaker + 1, order,
             weight);
        x->Decoder->SetOrderWeight(speaker, order, weight);
    } else if (method == "coeff") {
        if (argc < 5) {
            pd_error(x, "[ambi~]: Channel and Coeff not provided");
            return;
        }
        unsigned speaker = atom_getint(argv + 2) - 1;
        unsigned order = atom_getint(argv + 3);
        float weight = atom_getfloat(argv + 4);
        if (speaker > x->Decoder->GetSpeakerCount() - 1) {
            pd_error(x,
                     "[ambi~]: There is just %d speakers in the current "
                     "configuration",
                     x->Decoder->GetSpeakerCount());
            return;
        } else if (speaker < 0) {
            pd_error(x, "[ambi~]: Speaker must be a positive number");
            return;
        }

        if (order < 1) {
            pd_error(x, "[ambi~]: Order must be from 1 to 7");
            return;
        }

        post("[ambi~]: Speaker %d Order %d Weight %f", speaker + 1, order,
             weight);

        if (speaker > x->Decoder->GetSpeakerCount()) {
            pd_error(x,
                     "[ambi~]: There is just %d speakers in the configuration",
                     x->Decoder->GetSpeakerCount());
            return;
        } else if (speaker < 0) {
            pd_error(x, "[ambi~]: Speaker must be a positive number");
        }

        post("[ambi~]: Speaker %d Order %d Weight %f", speaker, order, weight);
        x->Decoder->SetCoefficient(speaker, order, weight);
    } else {
        pd_error(x, "[ambi~]: invalid method, options are: orderweight, coeff");
    }

    x->Decoder->Refresh();

    return;
}

// ─────────────────────────────────────
static void ambi_speakersconfig(ambi_tilde *x, t_symbol *s, int argc,
                                t_atom *argv) {
    std::string method = atom_getsymbol(argv + 1)->s_name;
    if (method == "pos") {
        if (argc < 3) {
            pd_error(x, "[ambi~]: HRTF file path not provided");
            return;
        }
        x->hrtfPath = atom_getsymbol(argv + 2)->s_name;
    }
    // set amount of speakers
    else if (method == "number") {
        if (!x->multichannel) {
            pd_error(x, "[ambi~]: Multichannel mode is off, you can't change "
                        "the amount of speakers");
            return;
        }
        if (argc < 3) {
            pd_error(x, "[ambi~]: Amount of speakers not provided");
            return;
        }
        unsigned nCh = atom_getint(argv + 2);
        if (nCh < 1) {
            pd_error(x,
                     "[ambi~]: Amount of speakers must be a positive number");
            return;
        }
        x->nChOut = nCh;
        post("[ambi~]: Amount of speakers set to %d", nCh);
        if (!ambi_speakerconfig(x)) {
            return;
        }
        canvas_update_dsp();
    }

    return;
}

// ─────────────────────────────────────
static void ambi_set(ambi_tilde *x, t_symbol *s, int argc, t_atom *argv) {
    if (argv[0].a_type != A_SYMBOL && argv[1].a_type != A_SYMBOL) {
        pd_error(x, "[ambi~] 1º and 2º must be a symbol");
        return;
    }

    std::string processor = atom_getsymbol(argv)->s_name;
    if (processor == "binaural") {
        ambi_binauralconfig(x, s, argc, argv);
    } else if (processor == "encoder") {
        ambi_encoderconfig(x, s, argc, argv);
    } else if (processor == "decoder") {
        ambi_decoderconfig(x, s, argc, argv);
    } else if (processor == "speakers") {
        ambi_speakersconfig(x, s, argc, argv);
    } else {
        pd_error(x, "[ambi~]: invalid processor, options are: binaural, "
                    "encoder, decoder");
    }
}

// ==============================================
static void ambi_azi(ambi_tilde *x, t_floatarg f) {
    x->Position.fAzimuth = DegreesToRadians(f);
}

// ─────────────────────────────────────
static void ambi_ele(ambi_tilde *x, t_floatarg f) {
    x->Position.fElevation = DegreesToRadians(f);
}

// ─────────────────────────────────────
static void ambi_dis(ambi_tilde *x, t_floatarg f) { x->Position.fDistance = f; }

// ─────────────────────────────────────
static t_int *ambi_perform(t_int *w) {

    ambi_tilde *x = (ambi_tilde *)(w[1]);
    unsigned n = (t_int)(w[2]);
    t_sample *in = (t_sample *)(w[3]);
    t_sample *out = (t_sample *)(w[4]);
    unsigned DspArr = 1 + 3 + x->nChOut;

    // silence the output
    if (!x->SpeakerPositionSet) {
        for (int i = 0; i < (x->nChOut * x->blockSize); i++) {
            out[i] = 0;
        }
        if (x->multichannel) {
            return w + 5;
        } else {
            return w + DspArr;
        }
    }

    unsigned int DecoderSpeakers = x->Decoder->GetSpeakerCount();
    x->Encoder->SetPosition(x->Position);
    x->Encoder->Process((float *)in, n, x->BFormat);

    if (x->binaural) {
        x->Binauralizer->Process(x->BFormat, (float **)x->SpeakersArray);
    } else {
        x->Decoder->Process(x->BFormat, n, (float **)x->SpeakersArray);
    }

    if (x->multichannel) {
        int outIndex = 0;
        for (int i = 0; i < x->nChOut; i++) {
            for (int j = 0; j < x->blockSize; j++) {
                out[outIndex] = x->SpeakersArray[i][j];
                outIndex++;
            }
        }
        printf("\n");
        printf("end copy\n");
        return (w + 5);
    } else {
        for (int i = 0; i < x->nChOut; i++) {
            t_sample *out = (t_sample *)(w[4 + i]);
            for (int j = 0; j < n; j++) {
                out[j] = x->SpeakersArray[i][j];
            }
        }
        return w + DspArr;
    }
}

// ─────────────────────────────────────
static void ambi_adddsp(ambi_tilde *x, t_signal **sp) {
    printf("ambi_adddsp\n");
    if (sp[0]->s_nchans != 1) {
        pd_error(x, "[ambi~]: ambi~ object just work with mono input!");
        return;
    }

    x->tailLength = 0;
    unsigned blockSize = sp[0]->s_n;
    unsigned sampleRate = sys_getsr();

    if (blockSize != x->blockSize || sampleRate != x->sampleRate ||
        x->Decoder->GetSpeakerCount() != x->nChOut) {
        ambi_configure(x, blockSize, sampleRate);
    }

    if (x->multichannel) {
        signal_setmultiout(&sp[1], x->nChOut);
        dsp_add(ambi_perform, 4, x, sp[0]->s_length * sp[0]->s_nchans,
                sp[0]->s_vec, sp[1]->s_vec);
    } else {
        for (int i = 0; i < x->nChOut; i++) {
            signal_setmultiout(&sp[i + 1], 1);
        }

        t_int *SigVec = (t_int *)getbytes((x->nChOut + 3) * sizeof(t_int));
        SigVec[0] = (t_int)x;
        SigVec[1] = (t_int)sp[0]->s_n;

        for (int i = 0; i < (x->nChOut + 1); i++) {
            SigVec[i + 2] = (t_int)sp[i]->s_vec;
        }
        dsp_addv(ambi_perform, x->nChOut + 3, SigVec);
        freebytes(SigVec, (x->nChOut + 3) * sizeof(t_int));
    }
    x->blockSize = blockSize;
}

// ─────────────────────────────────────
static void *ambi_new(t_symbol *s, int argc, t_atom *argv) {
    ambi_tilde *x = (ambi_tilde *)pd_new(ambi_class);
    x->xCanvas = canvas_getcurrent();
    t_symbol *patchDir = canvas_getdir(x->xCanvas);

    if (argv[0].a_type != A_FLOAT) {
        pd_error(nullptr,
                 "[ambi~]: First two arguments must be the number of input and "
                 "output channels");
        return NULL;
    }

    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            if (argv[i].a_type == A_SYMBOL) {
                std::string std = atom_getsymbol(argv + i)->s_name;
                if (std == "-m") {
                    x->multichannel = true;
                }
                if (std == "-b") {
                    x->binaural = true;
                }
            }
        }
    }

    x->exterDir = ambi_class->c_externdir->s_name;
    x->patchDir = patchDir->s_name;
    x->hrtfPath = x->exterDir + "/nh906.sofa";

    // Check if the file exists
    if (!fs::exists(x->hrtfPath)) {
        pd_error(x, "[ambi~]: Default HRTF not file found: %s",
                 x->hrtfPath.c_str());
    }

    x->nChOut = atom_getfloat(argv);
    x->outlets = (t_outlet **)malloc(x->nChOut * sizeof(t_outlet *));
    if (x->multichannel) {
        x->outlets[0] = outlet_new(&x->xObj, &s_signal);
    } else {
        for (int i = 0; i < x->nChOut; i++) {
            x->outlets[i] = outlet_new(&x->xObj, &s_signal);
        }
    }

    if (!ambi_speakerconfig(x)) {
        return nullptr;
    }

    if (x->NeedSpeakersPosition) {
        x->SpeakersSetup = new t_sample *[x->nChOut];
        for (int i = 0; i < x->nChOut; i++) {
            x->SpeakersSetup[i] = new t_sample[2];
            x->SpeakersSetup[i][0] = i;
            x->SpeakersSetup[i][1] = 0;
        }
    }

    // check how to define sample rate, needed?
    x->BFormat = new CBFormat();
    x->Encoder = new CAmbisonicEncoderDist();
    x->Decoder = new CAmbisonicDecoder();
    x->Binauralizer = new CAmbisonicBinauralizer();

    return x;
}

// ==============================================
static void *ambi_free(ambi_tilde *x) {
    delete x->BFormat;
    delete x->Encoder;
    delete x->Decoder;
    delete x->Binauralizer;

    unsigned int DecoderSpeakers = x->Decoder->GetSpeakerCount();

    for (int niSpeaker = 0; niSpeaker < DecoderSpeakers; niSpeaker++) {
        delete[] x->SpeakersArray[niSpeaker];
    }
    delete[] x->SpeakersArray;

    if (x->NeedSpeakersPosition) {
        for (int i = 0; i < x->nChOut; i++) {
            delete[] x->SpeakersSetup[i];
        }
    }

    free(x->outlets);
    return x;
}

// ─────────────────────────────────────
extern "C" void ambi_tilde_setup(void) {
    ambi_class = class_new(gensym("ambi~"), (t_newmethod)ambi_new,
                           (t_method)ambi_free, sizeof(ambi_tilde),
                           CLASS_DEFAULT | CLASS_MULTICHANNEL, A_GIMME, 0);

    CLASS_MAINSIGNALIN(ambi_class, ambi_tilde, sample);
    class_addmethod(ambi_class, (t_method)ambi_adddsp, gensym("dsp"), A_CANT,
                    0);

    // speaker, chNumber, Azi, Ele, Dis
    class_addmethod(ambi_class, (t_method)ambi_set, gensym("set"), A_GIMME, 0);

    // source position
    class_addmethod(ambi_class, (t_method)ambi_azi, gensym("azi"), A_FLOAT, 0);
    class_addmethod(ambi_class, (t_method)ambi_ele, gensym("ele"), A_FLOAT, 0);
    class_addmethod(ambi_class, (t_method)ambi_dis, gensym("dis"), A_FLOAT, 0);
}

#include <m_pd.h>

//
#include <string>

// libspatialaudio
#include <AmbisonicBinauralizer.h>
#include <AmbisonicDecoder.h>
#include <AmbisonicEncoderDist.h>
#include <BFormat.h>

static t_class *ambi;

// ==============================================
typedef struct _ambi {
    t_object xObj;
    t_canvas *xCanvas;
    t_sample sample;
    unsigned blockSize;
    unsigned sampleRate;

    unsigned nChIn;
    unsigned nChOut;
    std::string SpeakerConfig;
    Amblib_SpeakerSetUps SpeakerSetUp;
    bool NeedSpeakersPosition = false;
    bool SpeakerPositionSet = false;

    PolarPoint Position;
    t_float azimuth;
    t_float elevation;
    t_float distance;

    // Encoder config
    CAmbisonicEncoderDist *Encoder;

    // Decoder config
    CAmbisonicDecoder *Decoder;

    // Decoder config
    CBFormat *BFormat;

    // Binaural config
    CAmbisonicBinauralizer *Binauralizer;
    std::string HrtfPath;
    bool lowcpu;

    //
    t_sample **ppfSpeakerFeeds;
    t_sample **speakersSetup;

    bool DecoderConfigured = false;
    bool Binaural = false;

    t_inlet *Inlet;
    t_outlet **outlets;

} ambi_tilde;

// ==============================================
static int ambi_speakerconfig(ambi_tilde *x, std::string config) {
    if (config == "stereo") {
        x->SpeakerSetUp = Amblib_SpeakerSetUps::kAmblib_Stereo;
        x->SpeakerPositionSet = true;
        return 0;
    } else if (config == "quad") {
        x->SpeakerSetUp = Amblib_SpeakerSetUps::kAmblib_Quad;
        x->SpeakerPositionSet = true;
        return 0;
    } else if (config == "5.1") {
        x->SpeakerSetUp = Amblib_SpeakerSetUps::kAmblib_51;
        x->SpeakerPositionSet = true;
        return 0;
    } else if (config == "7.1") {
        x->SpeakerSetUp = Amblib_SpeakerSetUps::kAmblib_71;
        x->SpeakerPositionSet = true;
        return 0;
    } else if (config == "custom") {
        x->SpeakerSetUp = Amblib_SpeakerSetUps::kAmblib_CustomSpeakerSetUp;
        x->NeedSpeakersPosition = true;
        x->SpeakerPositionSet = false;
        return 0;
    } else {
        pd_error(x, "[ambi~]: invalid speaker configuration");
        x->SpeakerPositionSet = false;
        return -1;
    }
}

// ─────────────────────────────────────
static void ambi_binauralconfig(ambi_tilde *x, t_symbol *s, int argc,
                                t_atom *argv) {
    std::string method = atom_getsymbol(argv + 1)->s_name;
    if (method == "hrth") {
        if (argc < 3) {
            pd_error(x, "[ambi~]: HRTF file path not provided");
            return;
        }
        x->HrtfPath = atom_getsymbol(argv + 2)->s_name;
    } else if (method == "lowcpu") {
        bool lowcpu = atom_getfloat(argv + 2);
        if (lowcpu) {
            x->lowcpu = true;
        } else {
            x->lowcpu = false;
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
        if (argc < 5) {
            pd_error(x, "[ambi~]: Order and weight not provided");
            return;
        }
        float order = atom_getfloat(argv + 2);
        float weight = atom_getfloat(argv + 3);
        x->Encoder->SetOrderWeight(order, weight);
    } else if (method == "coeff") {
        if (argc < 5) {
            pd_error(x, "[ambi~]: Channel and Coeff not provided");
            return;
        }
        float channel = atom_getfloat(argv + 2);
        float coeff = atom_getfloat(argv + 3);
        x->Encoder->SetCoefficient(channel, coeff);
    } else if (method == "roomsize") {
        float size = atom_getfloat(argv + 2);
        x->Encoder->SetRoomRadius(size);
    }

    x->Encoder->Refresh();
}

// ─────────────────────────────────────
static void ambi_decoderconfig(ambi_tilde *x, t_symbol *s, int argc,
                               t_atom *argv) {
    std::string method = atom_getsymbol(argv + 1)->s_name;

    if (method == "orderweight") {
        if (argc < 5) {
            pd_error(x, "[ambi~]: Order and weight not provided");
            return;
        }

        float speaker = atom_getfloat(argv + 2);
        float order = atom_getfloat(argv + 3);
        float weight = atom_getfloat(argv + 4);
        x->Decoder->SetOrderWeight(speaker, order, weight);
    } else if (method == "coeff") {
        if (argc < 5) {
            pd_error(x, "[ambi~]: Channel and Coeff not provided");
            return;
        }
        float speaker = atom_getfloat(argv + 2);
        float order = atom_getfloat(argv + 3);
        float weight = atom_getfloat(argv + 4);
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
        x->HrtfPath = atom_getsymbol(argv + 2)->s_name;
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
        ambi_encoderconfig(x, s, argc, argv);
    } else if (processor == "speakers") {
        ambi_speakersconfig(x, s, argc, argv);
    } else {
        pd_error(x, "[ambi~]: invalid processor, options are: binaural, "
                    "encoder, decoder");
    }
}

// ==============================================
static void ambi_speakers(ambi_tilde *x, t_int nCh, t_float azi, t_float ele,
                          t_float dis) {

    nCh = (int)nCh - 1; // fix index for the channels numbers

    if (!x->DecoderConfigured) {
        pd_error(x, "[ambi~]: Decoder not configured yet, turn the DSP on");
        return;
    }

    if (nCh > x->Decoder->GetSpeakerCount()) {
        pd_error(x, "[ambi~]: There is just %d speakers in the configuration",
                 x->Decoder->GetSpeakerCount());
    }

    PolarPoint newSpeakerPosition;
    newSpeakerPosition.fAzimuth = DegreesToRadians(azi);
    newSpeakerPosition.fElevation = DegreesToRadians(ele);
    newSpeakerPosition.fDistance = dis;
    x->Decoder->SetPosition(nCh, newSpeakerPosition);

    x->speakersSetup[nCh][1] = 1;
    post("[ambi~]: Speaker channel %d set", (int)nCh + 1);

    bool allSpeakersSet = true;
    for (int i = 0; i < x->nChOut; i++) {
        if (x->speakersSetup[i][1] == 0) {
            allSpeakersSet = false;
            break;
        }
    }
    if (allSpeakersSet && !x->SpeakerPositionSet) {
        x->SpeakerPositionSet = true;
        post("[ambi~]: All speakers set");
    } else {
        x->Decoder->Refresh();
    }
    return;
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
    unsigned n = (int)(w[2]);
    unsigned DspArr = 1 + x->nChOut + 3;
    unsigned outChIndex = 3 + 1;

    for (int i = 0; i < x->nChOut; i++) {
        post("float %f", w[i]);
    }

    // silence the output
    if (!x->SpeakerPositionSet) {
        for (int i = 0; i < x->nChOut; i++) {
            t_sample *out = (t_sample *)(w[outChIndex + i]);
            for (int j = 0; j < n; j++) {
                out[j] = 0;
            }
        }
        return w + DspArr;
    }

    // Position
    unsigned int DecoderSpeakers = x->Decoder->GetSpeakerCount();
    x->Encoder->SetPosition(x->Position);

    // Encoder
    t_sample *in = (t_sample *)(w[3]);
    x->Encoder->Process(in, n, x->BFormat);

    // Process the audio
    if (x->Binaural) {
        x->Binauralizer->Process(x->BFormat, x->ppfSpeakerFeeds);
    } else {
        x->Decoder->Process(x->BFormat, n, x->ppfSpeakerFeeds);
    }

    // copy output to the out
    for (int i = 0; i < x->nChOut; i++) {
        t_sample *out = (t_sample *)(w[outChIndex + i]);
        for (int j = 0; j < n; j++) {
            out[j] = x->ppfSpeakerFeeds[i][j];
        }
    }

    return w + DspArr;
}

// ==============================================
static void ambi_adddsp(ambi_tilde *x, t_signal **sp) {
    unsigned tailLength;
    int unsigned blockSize = sp[0]->s_n;
    int unsigned sampleRate = sys_getsr();

    if (blockSize != x->blockSize || sampleRate != x->sampleRate) {
        x->BFormat->Configure(1, true, blockSize);

        float fadeTimeInMilliSec =
            1000.f * (float)blockSize / (float)sampleRate;

        x->Encoder->Configure(1, true, sys_getsr());
        x->Decoder->Configure(1, true, blockSize, sys_getsr(), x->SpeakerSetUp,
                              x->nChOut);
        x->Binauralizer->Configure(1, true, sys_getsr(), blockSize, tailLength,
                                   x->HrtfPath);

        // memory for speaker feeds
        unsigned int DecoderSpeakers = x->Decoder->GetSpeakerCount();
        x->ppfSpeakerFeeds = new float *[DecoderSpeakers];
        for (int niSpeaker = 0; niSpeaker < DecoderSpeakers; niSpeaker++) {
            x->ppfSpeakerFeeds[niSpeaker] = new float[blockSize];
        }

        x->blockSize = blockSize;
        x->sampleRate = sampleRate;
        x->DecoderConfigured = true;
    }

    // this is from circuit~ :)
    t_int *SigVec = (t_int *)getbytes((x->nChOut + 3) * sizeof(t_int));
    SigVec[0] = (t_int)x;
    SigVec[1] = (t_int)sp[0]->s_n;

    for (int i = 0; i < (x->nChOut + 1); i++) {
        SigVec[i + 2] = (t_int)sp[i]->s_vec;
    }

    dsp_addv(ambi_perform, x->nChOut + 3, SigVec);
    freebytes(SigVec, (x->nChOut + 3) * sizeof(t_int));
}

// ─────────────────────────────────────
static void *ambi_new(t_symbol *s, int argc, t_atom *argv) {
    ambi_tilde *x = (ambi_tilde *)pd_new(ambi);
    x->xCanvas = canvas_getcurrent();
    t_symbol *patchDir = canvas_getdir(x->xCanvas);

    // check if argv[0] and argc[1] are t_float
    // first two are n input chns and output
    if (argv[0].a_type != A_FLOAT) {
        pd_error(nullptr,
                 "[ambi~]: First two arguments must be the number of input and "
                 "output channels");
        return NULL;
    }
    x->nChOut = atom_getfloat(argv);

    if (argc > 2) {
        for (int i = 2; i < argc; i++) {
            if (argv[i].a_type == A_SYMBOL) {
                std::string thing = atom_getsymbol(argv + i)->s_name;
                if (thing == "-s") {
                    i++;
                    if (argv[i].a_type == A_SYMBOL) {
                        x->SpeakerConfig = atom_getsymbol(argv + i)->s_name;
                    } else if (argv[i].a_type == A_FLOAT) {
                        x->SpeakerConfig =
                            std::to_string(atom_getfloat(argv + i));
                    }
                } else if (thing == "-b") {
                    x->Binaural = true;
                    x->nChOut = 2;
                }
            }
        }
    }
    if (x->Binaural && x->HrtfPath.empty()) {
        post("[ambi~] Using default HRTF file");
    }

    // malloc outlets for nChOut
    x->outlets = (t_outlet **)malloc(x->nChOut * sizeof(t_outlet *));
    for (int i = 0; i < x->nChOut; i++) {
        x->outlets[i] = outlet_new(&x->xObj, &s_signal);
    }

    if (x->SpeakerConfig.empty()) {
        if (x->nChOut == 2) {
            x->SpeakerConfig = "stereo";
        } else if (x->nChOut == 4) {
            x->SpeakerConfig = "quad";
        } else {
            pd_error(x, "[ambi~]: You are using a custom speaker configuration "
                        "but you haven't specified the number of speakers. Use "
                        "[-s custom] to set.");
            return NULL;
        }
    }

    ambi_speakerconfig(x, x->SpeakerConfig);
    if (x->NeedSpeakersPosition) {
        // create an array with the channel number and if it was configured or
        // not
        x->speakersSetup = new float *[x->nChOut];
        for (int i = 0; i < x->nChOut; i++) {
            x->speakersSetup[i] = new float[2];
            x->speakersSetup[i][0] = i;
            x->speakersSetup[i][1] = 0;
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
        delete[] x->ppfSpeakerFeeds[niSpeaker];
    }
    delete[] x->ppfSpeakerFeeds;

    if (x->NeedSpeakersPosition) {
        for (int i = 0; i < x->nChOut; i++) {
            delete[] x->speakersSetup[i];
        }
    }

    free(x->outlets);
    return (void *)x;
}

// ─────────────────────────────────────
extern "C" void ambi_tilde_setup(void) {
    ambi =
        class_new(gensym("ambi~"), (t_newmethod)ambi_new, (t_method)ambi_free,
                  sizeof(ambi_tilde), CLASS_DEFAULT, A_GIMME, 0);

    CLASS_MAINSIGNALIN(ambi, ambi_tilde, sample);
    class_addmethod(ambi, (t_method)ambi_adddsp, gensym("dsp"), A_CANT, 0);

    // speaker, chNumber, Azi, Ele, Dis
    class_addmethod(ambi, (t_method)ambi_speakers, gensym("speaker"), A_FLOAT,
                    A_FLOAT, A_FLOAT, A_FLOAT, 0);
    class_addmethod(ambi, (t_method)ambi_set, gensym("set"), A_GIMME, 0);

    // source position
    class_addmethod(ambi, (t_method)ambi_azi, gensym("azi"), A_FLOAT, 0);
    class_addmethod(ambi, (t_method)ambi_ele, gensym("ele"), A_FLOAT, 0);
    class_addmethod(ambi, (t_method)ambi_dis, gensym("dis"), A_FLOAT, 0);
}

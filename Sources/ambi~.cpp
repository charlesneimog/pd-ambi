#include <m_pd.h>

//
#include <fstream>
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
    CBFormat *BFormat;

    CAmbisonicEncoderDist *Encoder;
    CAmbisonicDecoder *Decoder;
    CAmbisonicBinauralizer *Binauralizer;

    std::string HrtfPath;
    float **ppfSpeakerFeeds;

    float **speakersSetup;

    bool DecoderConfigured = false;
    bool Binaural = false;

    t_float azimuth;
    t_float elevation;
    t_float distance;

    t_inlet **inlets;
    t_outlet **outlets;

} elseAmbi;

// ==============================================
static int SetSpeakerConfig(elseAmbi *x, std::string config) {
    if (config == "stereo") {
        x->SpeakerSetUp = kAmblib_Stereo;
        x->SpeakerPositionSet = true;
        return 0;
    } else if (config == "quad") {
        x->SpeakerSetUp = kAmblib_Quad;
        x->SpeakerPositionSet = true;
        return 0;
    } else if (config == "5.1") {
        x->SpeakerSetUp = kAmblib_51;
        x->SpeakerPositionSet = true;
        return 0;
    } else if (config == "7.1") {
        x->SpeakerSetUp = kAmblib_71;
        x->SpeakerPositionSet = true;
        return 0;
    } else if (config == "custom") {
        x->SpeakerSetUp = kAmblib_CustomSpeakerSetUp;
        x->NeedSpeakersPosition = true;
        x->SpeakerPositionSet = false;
        return 0;
    } else {
        pd_error(x, "[ambi~]: invalid speaker configuration");
        x->SpeakerPositionSet = false;
        return -1;
    }
}

// ==============================================
static void NewSpeaker(elseAmbi *x, t_float nCh, t_float azi, t_float ele,
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

    x->speakersSetup[(int)nCh][1] = 1;
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
static void SetAzimuth(elseAmbi *x, t_floatarg f) {
    x->Position.fAzimuth = DegreesToRadians(f);
}
static void SetElevation(elseAmbi *x, t_floatarg f) {
    x->Position.fElevation = DegreesToRadians(f);
}
static void SetDistance(elseAmbi *x, t_floatarg f) {
    x->Position.fDistance = f;
}

// ==============================================
static void ToggleBinaural(elseAmbi *x, t_floatarg f) {
    if (f == 0) {
        x->Binaural = false;
        post("[ambi~]: Binaural mode disabled");
    } else {
        x->Binaural = true;
        post("[ambi~]: Binaural mode enabled");
    }
}

// ==============================================
static t_int *AmbiPerform(t_int *w) {
    elseAmbi *x = (elseAmbi *)(w[1]);
    unsigned n = (int)(w[2]);
    unsigned DspArr = x->nChIn + x->nChOut + 3;
    unsigned inChIndex = 3;
    unsigned outChIndex = 3 + x->nChIn;

    if (!x->SpeakerPositionSet) {
        for (int i = 0; i < x->nChOut; i++) {
            t_sample *out = (t_sample *)(w[outChIndex + i]);
            for (int j = 0; j < n; j++) {
                out[j] = 0;
            }
        }
    }

    unsigned int DecoderSpeakers = x->Decoder->GetSpeakerCount();
    x->Encoder->SetPosition(x->Position); // source position

    for (int i = 0; i < x->nChIn; i++) {
        t_sample *in = (t_sample *)(w[inChIndex + i]);
        x->Encoder->Process(in, n, x->BFormat);
    }

    if (x->Binaural) {
        x->Binauralizer->Process(x->BFormat, x->ppfSpeakerFeeds);
    } else {
        x->Decoder->Process(x->BFormat, n, x->ppfSpeakerFeeds);
    }
    for (int i = 0; i < x->nChOut; i++) {
        t_sample *out = (t_sample *)(w[outChIndex + i]);
        for (int j = 0; j < n; j++) {
            out[j] = x->ppfSpeakerFeeds[i][j];
        }
    }

    return w + DspArr;
}

// ==============================================
static void AmbiAddDsp(elseAmbi *x, t_signal **sp) {

    unsigned int ChCount = x->nChIn + x->nChOut;
    unsigned tailLength;
    int unsigned blockSize = sp[0]->s_n;
    int unsigned sampleRate = sys_getsr();

    // this must be done once or when something change
    if (blockSize != x->blockSize || sampleRate != x->sampleRate) {
        x->BFormat->Configure(1, true, blockSize); // set
        x->Encoder->Configure(1, true, sys_getsr());
        x->Decoder->Configure(1, true, blockSize, x->SpeakerSetUp, x->nChOut);
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
    t_int *SigVec = (t_int *)getbytes((ChCount + 2) * sizeof(t_int));
    SigVec[0] = (t_int)x;
    SigVec[1] = (t_int)sp[0]->s_n;

    for (int i = 0; i < ChCount; i++) {
        SigVec[i + 2] = (t_int)sp[i]->s_vec;
    }

    dsp_addv(AmbiPerform, ChCount + 2, SigVec);
    freebytes(SigVec, (ChCount + 2) * sizeof(t_int));
}

// ─────────────────────────────────────
static void *NewAmbi(t_symbol *s, int argc, t_atom *argv) {
    elseAmbi *x = (elseAmbi *)pd_new(ambi);
    x->xCanvas = canvas_getcurrent();
    t_symbol *patchDir = canvas_getdir(x->xCanvas);

    // check if argv[0] and argc[1] are t_float
    // first two are n input chns and output
    if (argv[0].a_type != A_FLOAT || argv[1].a_type != A_FLOAT) {
        pd_error(nullptr,
                 "[ambi~]: First two arguments must be the number of input and "
                 "output channels");
        return NULL;
    }

    if (atom_getfloat(argv) != 1) {
        pd_error(x, "[ambi~]: Input channels different from 1 are not "
                    "supported yet!");
        return NULL;
    }

    x->nChIn = atom_getfloat(argv);
    x->nChOut = atom_getfloat(argv + 1);

    // malloc inlets for nChIn
    x->inlets = (t_inlet **)malloc(x->nChIn * sizeof(t_inlet *));
    for (int i = 0; i < x->nChIn - 1; i++) {
        x->inlets[i] =
            inlet_new(&x->xObj, &x->xObj.ob_pd, &s_signal, &s_signal);
    }

    if (argc > 2) {
        for (int i = 2; i < argc; i++) {
            if (argv[i].a_type == A_SYMBOL) {
                std::string thing = atom_getsymbol(argv + i)->s_name;
                if (thing == "-s") { // find better name
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
                } else if (thing == "-hrtf") {
                    if (i + 1 < argc) {
                        i++;
                        if (argv[i].a_type == A_SYMBOL) {
                            std::string hrtfPath =
                                atom_getsymbol(argv + i)->s_name;
                            if (hrtfPath.size() < 5) {
                                hrtfPath += ".sofa";
                            } else {
                                if (hrtfPath.substr(hrtfPath.size() - 5) ==
                                    ".sofa") {
                                    x->HrtfPath = hrtfPath;
                                } else {
                                    hrtfPath += ".sofa";
                                }
                            }
                            // patchDir + hrtfPath
                            std::string completePath = patchDir->s_name;
                            completePath += "/";
                            completePath += hrtfPath;
                            std::ifstream file(completePath);
                            if (file.good()) {
                                x->HrtfPath = completePath;
                            } else {
                                pd_error(x, "[ambi~]: HRTF file not found, "
                                            "using default");
                            }
                        }
                    }
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

    SetSpeakerConfig(x, x->SpeakerConfig);
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
static void *FreeAmbi(elseAmbi *x) {
    delete x->BFormat;
    delete x->Encoder;
    delete x->Decoder;
    delete x->Binauralizer;

    unsigned int DecoderSpeakers = x->Decoder->GetSpeakerCount();

    for (int niSpeaker = 0; niSpeaker < DecoderSpeakers; niSpeaker++) {
        delete[] x->ppfSpeakerFeeds[niSpeaker];
    }
    delete[] x->ppfSpeakerFeeds;

    free(x->outlets);
    return (void *)x;
}

extern "C" void ambi_tilde_setup(void);

// ─────────────────────────────────────
void ambi_tilde_setup(void) {
    ambi = class_new(gensym("ambi~"), (t_newmethod)NewAmbi, (t_method)FreeAmbi,
                     sizeof(elseAmbi), CLASS_DEFAULT, A_GIMME,
                     0); // maybe change to A_FLOAT AFLOAT, A_GIME

    CLASS_MAINSIGNALIN(ambi, elseAmbi, sample);
    class_addmethod(ambi, (t_method)AmbiAddDsp, gensym("dsp"), A_CANT, 0);

    // speaker, chNumber, Azi, Ele, Dis
    class_addmethod(ambi, (t_method)NewSpeaker, gensym("speaker"), A_FLOAT,
                    A_FLOAT, A_FLOAT, A_FLOAT, 0);

    // source position
    class_addmethod(ambi, (t_method)SetAzimuth, gensym("azi"), A_FLOAT, 0);
    class_addmethod(ambi, (t_method)SetElevation, gensym("ele"), A_FLOAT, 0);
    class_addmethod(ambi, (t_method)SetDistance, gensym("dis"), A_FLOAT, 0);

    // another things
    class_addmethod(ambi, (t_method)ToggleBinaural, gensym("binaural"), A_FLOAT,
                    0);
}

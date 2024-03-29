#include <m_pd.h>

//
#include <string>

// libspatialaudio
#include <AmbisonicDecoder.h>
#include <AmbisonicEncoder.h>
#include <BFormat.h>

static t_class *ambi;

// ─────────────────────────────────────
typedef struct _ambi {
    t_object xObj;
    t_sample sample;

    unsigned int nChIn;
    unsigned int nChOut;
    std::string SpeakerConfig;
    Amblib_SpeakerSetUps SpeakerSetUp;
    bool NeedSpeakersPosition = false;
    bool SpeakerPositionSet = false;

    PolarPoint Position;
    CBFormat *BFormat;
    CAmbisonicEncoder *Encoder;
    CAmbisonicDecoder *Decoder;

    bool DecoderConfigured = false;
    bool Binaural = true;

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

    if (nCh != x->nChOut) {
        pd_error(x, "[ambi~]: The number of channels must be equal to the "
                    "number of speakers in the configuration");
        return;
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
static t_int *AmbiPerform(t_int *w) {
    elseAmbi *x = (elseAmbi *)(w[1]);
    unsigned int n = (int)(w[2]);

    unsigned int DspArr = x->nChIn + x->nChOut + 3;

    unsigned int inChIndex = 3;
    unsigned int outChIndex = 3 + x->nChIn;

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

    float **ppfSpeakerFeeds = new float *[DecoderSpeakers];
    for (int niSpeaker = 0; niSpeaker < DecoderSpeakers; niSpeaker++) {
        ppfSpeakerFeeds[niSpeaker] = new float[n];
    }
    x->Decoder->Process(x->BFormat, n, ppfSpeakerFeeds);

    for (int i = 0; i < x->nChOut; i++) {
        t_sample *out = (t_sample *)(w[outChIndex + i]);
        for (int j = 0; j < n; j++) {
            out[j] = ppfSpeakerFeeds[i][j];
        }
    }

    for (int niSpeaker = 0; niSpeaker < DecoderSpeakers; niSpeaker++) {
        delete[] ppfSpeakerFeeds[niSpeaker];
    }
    delete[] ppfSpeakerFeeds;

    return w + DspArr;
}

// ==============================================
static void AmbiAddDsp(elseAmbi *x, t_signal **sp) {

    unsigned int ChCount = x->nChIn + x->nChOut;

    int unsigned blockSize = sp[0]->s_n;
    x->BFormat->Configure(1, true, blockSize);
    x->Encoder->Configure(1, true, 0);
    x->Decoder->Configure(1, true, blockSize, x->SpeakerSetUp, x->nChOut);

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

    // check if argv[0] and argc[1] are t_float
    // first two are n input chns and output
    if (argv[0].a_type != A_FLOAT || argv[1].a_type != A_FLOAT) {
        post("[ambi~]: invalid arguments");
        return NULL;
    }

    if (atom_getfloat(argv) != 1) {
        pd_error(x, "[ambi~]: Input channels different from 1 are not "
                    "supported yet");
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

    if (argc > 3) {
        for (int i = 3; i < argc; i++) {
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
                }
            }
        }
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

    // check how to define sample rate, needed?
    x->BFormat = new CBFormat();
    x->Encoder = new CAmbisonicEncoder();
    x->Decoder = new CAmbisonicDecoder();

    return x;
}

static void *FreeAmbi(elseAmbi *x) {
    delete x->BFormat;
    delete x->Encoder;
    delete x->Decoder;
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

    class_addmethod(ambi, (t_method)SetAzimuth, gensym("azi"), A_FLOAT, 0);
    class_addmethod(ambi, (t_method)SetElevation, gensym("ele"), A_FLOAT, 0);
    class_addmethod(ambi, (t_method)SetDistance, gensym("dis"), A_FLOAT, 0);
}

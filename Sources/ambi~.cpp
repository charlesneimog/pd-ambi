#include <m_pd.h>

// libspatialaudio
#include <AmbisonicDecoder.h>
#include <AmbisonicEncoder.h>
#include <BFormat.h>

static t_class *ambi;

// ─────────────────────────────────────
typedef struct _ambi {
    t_object xObj;
    t_sample sample;

    PolarPoint Position;
    CBFormat *BFormat;
    CAmbisonicEncoder *Encoder;
    CAmbisonicDecoder *Decoder;

    t_float azimuth;
    t_float elevation;
    t_float distance;

    t_outlet *out1;
    t_outlet *out2;
    t_outlet *out3;
    t_outlet *out4;

} elseAmbi;

static void SetAzimuth(elseAmbi *x, t_floatarg f) { x->Position.fAzimuth = f; }
static void SetElevation(elseAmbi *x, t_floatarg f) {
    x->Position.fElevation = f;
}
static void SetDistance(elseAmbi *x, t_floatarg f) {
    x->Position.fDistance = f;
}

// ==============================================
static t_int *AmbiPerform(t_int *w) {
    elseAmbi *x = (elseAmbi *)(w[1]);
    t_sample *in = (t_sample *)(w[2]);

    t_sample *out1 = (t_sample *)(w[3]);
    t_sample *out2 = (t_sample *)(w[4]);
    t_sample *out3 = (t_sample *)(w[5]);
    t_sample *out4 = (t_sample *)(w[6]);

    int n = (int)(w[7]);

    unsigned int DecoderSpeakers = x->Decoder->GetSpeakerCount();
    x->Encoder->SetPosition(x->Position);

    x->Encoder->Process(in, n, x->BFormat);

    // rethink
    float **ppfSpeakerFeeds = new float *[DecoderSpeakers];
    for (int niSpeaker = 0; niSpeaker < DecoderSpeakers; niSpeaker++) {
        ppfSpeakerFeeds[niSpeaker] = new float[n];
    }

    x->Decoder->Process(x->BFormat, n, ppfSpeakerFeeds);

    // just stereo for test
    for (int i = 0; i < n; i++) {
        out1[i] = ppfSpeakerFeeds[0][i];
        out2[i] = ppfSpeakerFeeds[1][i];
        out3[i] = ppfSpeakerFeeds[2][i];
        out4[i] = ppfSpeakerFeeds[3][i];
    }

    for (int niSpeaker = 0; niSpeaker < DecoderSpeakers; niSpeaker++)
        delete[] ppfSpeakerFeeds[niSpeaker];
    delete[] ppfSpeakerFeeds;

    return (w + 8);
}

// ==============================================
static void AmbiAddDsp(elseAmbi *x, t_signal **sp) {
    // bool CAmbisonicDecoder::Configure(unsigned nOrder, bool b3D,
    //                                   unsigned nBlockSize, int nSpeakerSetUp,
    //                                   unsigned nSpeakers) {
    //

    int unsigned blockSize = sp[0]->s_n;
    x->BFormat->Configure(1, true, blockSize);
    x->Encoder->Configure(1, true, 0); // nao sei oq é esse terceiro parametro
    x->Decoder->Configure(1, true, blockSize, kAmblib_Cube,
                          4); // 4 é o numero de canais
    x->Encoder->SetPosition(x->Position);

    // tem que variar dependendo do que o usuario vai querer
    dsp_add(AmbiPerform, 7, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec,
            sp[3]->s_vec, sp[4]->s_vec, (t_int)sp[0]->s_n);
}

// ─────────────────────────────────────
static void *NewAmbi(t_symbol *s, int argc, t_atom *argv) {
    elseAmbi *x = (elseAmbi *)pd_new(ambi);

    x->out1 = outlet_new(&x->xObj, &s_signal);
    x->out2 = outlet_new(&x->xObj, &s_signal);
    x->out3 = outlet_new(&x->xObj, &s_signal);
    x->out4 = outlet_new(&x->xObj, &s_signal);

    // check how to define sample rate, needed?
    x->BFormat = new CBFormat();
    x->Encoder = new CAmbisonicEncoder();
    x->Decoder = new CAmbisonicDecoder();

    return x;
}

extern "C" void ambi_tilde_setup(void);

// ─────────────────────────────────────
void ambi_tilde_setup(void) {
    ambi = class_new(gensym("ambi~"), (t_newmethod)NewAmbi, NULL,
                     sizeof(elseAmbi), CLASS_DEFAULT, A_GIMME, 0);

    CLASS_MAINSIGNALIN(ambi, elseAmbi, sample);
    class_addmethod(ambi, (t_method)AmbiAddDsp, gensym("dsp"), A_CANT, 0);

    class_addmethod(ambi, (t_method)SetAzimuth, gensym("azi"), A_FLOAT, 0);
    class_addmethod(ambi, (t_method)SetElevation, gensym("ele"), A_FLOAT, 0);
    class_addmethod(ambi, (t_method)SetDistance, gensym("dis"), A_FLOAT, 0);
}

#include <assert.h>
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "api.h"
#include "randombytes.h"

#if MASKED
#include "masked_api.h"
#endif

#if REQUIRES_BUF
#define BUF_BYTES 1000000
#endif

// https://stackoverflow.com/a/1489985/1711232
#define PASTER(x, y) x##_##y
#define EVALUATOR(x, y) PASTER(x, y)
#define NAMESPACE(fun) EVALUATOR(PQCLEAN_NAMESPACE, fun)

#define CRYPTO_PUBLICKEYBYTES NAMESPACE(CRYPTO_PUBLICKEYBYTES)
#define CRYPTO_SECRETKEYBYTES NAMESPACE(CRYPTO_SECRETKEYBYTES)
#define CRYPTO_BYTES          NAMESPACE(CRYPTO_BYTES)
#define CRYPTO_ALGNAME        NAMESPACE(CRYPTO_ALGNAME)

#define crypto_sign_keypair NAMESPACE(crypto_sign_keypair)
#define crypto_sign NAMESPACE(crypto_sign)
#define crypto_sign_open NAMESPACE(crypto_sign_open)

#if REQUIRES_BUF
#define GEN_BUF_BYTES NAMESPACE(GEN_BUF_BYTES)
#define SIGN_BUF_BYTES NAMESPACE(SIGN_BUF_BYTES)
#define VER_BUF_BYTES NAMESPACE(VER_BUF_BYTES)
#endif

#if MASKED
#define unmask_sk NAMESPACE(unmask_sk)
#endif

void nist_kat_init(unsigned char *entropy_input, unsigned char *personalization_string, int security_strength);

static void fprintBstr(FILE *fp, const char *S, const uint8_t *A, size_t L) {
    size_t i;
    fprintf(fp, "%s", S);
    for (i = 0; i < L; i++) {
        fprintf(fp, "%02X", A[i]);
    }
    if (L == 0) {
        fprintf(fp, "00");
    }
    fprintf(fp, "\n");
}

int main() {

    uint8_t entropy_input[48];
    uint8_t seed[48];
    FILE *fh = stdout;
    uint8_t public_key[CRYPTO_PUBLICKEYBYTES];
    uint8_t secret_key[CRYPTO_SECRETKEYBYTES];
    size_t mlen = 33;
    size_t smlen, mlen1;
    uint8_t m[33];
    uint8_t sm[33 + CRYPTO_BYTES];
    int rc;
#if REQUIRES_BUF
    uint8_t gen_buf[GEN_BUF_BYTES];
    uint8_t sig_buf[SIGN_BUF_BYTES];
    uint8_t ver_buf[VER_BUF_BYTES];
#endif

#if MASKED
	uint8_t umsk_secret_key[CRYPTO_SECRETKEYBYTES];
	int umsk_crypto_secretbytes = 0;
#endif
    for (uint8_t i = 0; i < 48; i++) {
        entropy_input[i] = i;
    }

    nist_kat_init(entropy_input, NULL, 256);

    fprintf(fh, "count = 0\n");
    randombytes(seed, 48);
    fprintBstr(fh, "seed = ", seed, 48);

    fprintf(fh, "mlen = 33\n");

    randombytes(m, mlen);
    fprintBstr(fh, "msg = ", m, mlen);

    nist_kat_init(seed, NULL, 256);

#if REQUIRES_BUF
    rc = crypto_sign_keypair(public_key, secret_key, gen_buf);
#else 
    rc = crypto_sign_keypair(public_key, secret_key);
#endif
    if (rc != 0) {
        fprintf(stderr, "[kat_kem] %s ERROR: crypto_kem_keypair failed!\n", CRYPTO_ALGNAME);
        return -1;
    }
    fprintBstr(fh, "pk = ", public_key, CRYPTO_PUBLICKEYBYTES);
#if MASKED
	umsk_crypto_secretbytes = unmask_sk(umsk_secret_key, secret_key);
    fprintBstr(fh, "sk = ", umsk_secret_key, umsk_crypto_secretbytes);
#else
    fprintBstr(fh, "sk = ", secret_key, CRYPTO_SECRETKEYBYTES);
#endif

#if REQUIRES_BUF
    rc = crypto_sign(sm, &smlen, m, mlen, secret_key, sig_buf);
#else
	rc = crypto_sign(sm, &smlen, m, mlen, secret_key);
#endif
    if (rc != 0) {
        fprintf(stderr, "[kat_kem] %s ERROR: crypto_sign failed!\n", CRYPTO_ALGNAME);
        return -2;
    }
    fprintf(fh, "smlen = %zu\n", smlen);
    fprintBstr(fh, "sm = ", sm, smlen);

#if REQUIRES_BUF
    rc = crypto_sign_open(sm, &mlen1, sm, smlen, public_key, ver_buf);
#else
	rc = crypto_sign_open(sm, &mlen1, sm, smlen, public_key);
#endif
    if (rc != 0) {
        fprintf(stderr, "[kat_kem] %s ERROR: crypto_sign_open failed!\n", CRYPTO_ALGNAME);
        return -3;
    }

    if ( mlen != mlen1 ) {
        printf("crypto_sign_open returned bad 'mlen': got <%zu>, expected <%zu>\n", mlen1, mlen);
        return -4;
    }
    if (memcmp(m, sm, mlen)) {
        printf("crypto_sign_open returned bad 'm' value\n");
        return -5;
    }
    return 0;

}

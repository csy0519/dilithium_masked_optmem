#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "api.h"
#include "randombytes.h"

#if MASKED
#include "masked_api.h"
#endif

#define MAXMLEN 2048

#if REQUIRES_BUF
#define BUF_BYTES 1000000
#endif


static void printbytes(const uint8_t *x, size_t xlen) {
    size_t i;
    for (i = 0; i < xlen; i++) {
        printf("%02x", x[i]);
    }
    printf("\n");
}

// https://stackoverflow.com/a/1489985/1711232
#define PASTER(x, y) x##_##y
#define EVALUATOR(x, y) PASTER(x, y)
#define NAMESPACE(fun) EVALUATOR(PQCLEAN_NAMESPACE, fun)

#define CRYPTO_PUBLICKEYBYTES NAMESPACE(CRYPTO_PUBLICKEYBYTES)
#define CRYPTO_SECRETKEYBYTES NAMESPACE(CRYPTO_SECRETKEYBYTES)
#define CRYPTO_BYTES          NAMESPACE(CRYPTO_BYTES)

#define crypto_sign_keypair NAMESPACE(crypto_sign_keypair)
#define crypto_sign NAMESPACE(crypto_sign)
#define crypto_sign_open NAMESPACE(crypto_sign_open)
#define crypto_sign_signature NAMESPACE(crypto_sign_signature)
#define crypto_sign_verify NAMESPACE(crypto_sign_verify)

#if REQUIRES_BUF
#define GEN_BUF_BYTES NAMESPACE(GEN_BUF_BYTES)
#define SIGN_BUF_BYTES NAMESPACE(SIGN_BUF_BYTES)
#define VER_BUF_BYTES NAMESPACE(VER_BUF_BYTES)
#endif

#if MASKED
#define unmask_sk NAMESPACE(unmask_sk)
#endif

int main(void) {
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];

    uint8_t mi[MAXMLEN];
    uint8_t sm[MAXMLEN + CRYPTO_BYTES];
    uint8_t sig[CRYPTO_BYTES];

#if REQUIRES_BUF
    uint8_t gen_buf[GEN_BUF_BYTES];
    uint8_t sig_buf[SIGN_BUF_BYTES];
    uint8_t ver_buf[VER_BUF_BYTES];
#endif

#if MASKED
	uint8_t umsk_sk[CRYPTO_SECRETKEYBYTES];
	int umsk_crypto_secretbytes = 0;
#endif

    size_t smlen;
    size_t siglen;
    size_t mlen;

    int r;
    size_t i, k;

    /* i = 0, 1, 4, 16, 64, 256, 1024 */
    for (i = 0; i < MAXMLEN; i = (i == 0) ? i + 1 : i << 2) {
        randombytes(mi, i);

#if REQUIRES_BUF
        crypto_sign_keypair(pk, sk, gen_buf);
#else 
        crypto_sign_keypair(pk, sk);
#endif
        printbytes(pk, CRYPTO_PUBLICKEYBYTES);

#if MASKED
        umsk_crypto_secretbytes = unmask_sk(umsk_sk, sk);
		printbytes(umsk_sk, umsk_crypto_secretbytes);
#else
		printbytes(sk, CRYPTO_SECRETKEYBYTES);
#endif

#if REQUIRES_BUF
        crypto_sign(sm, &smlen, mi, i, sk, sig_buf);
        crypto_sign_signature(sig, &siglen, mi, i, sk, sig_buf);
#else 
        crypto_sign(sm, &smlen, mi, i, sk);
        crypto_sign_signature(sig, &siglen, mi, i, sk);
#endif
        printbytes(sm, smlen);
        printbytes(sig, siglen);

        // By relying on m == sm we prevent having to allocate CRYPTO_BYTES
        // twice
#if REQUIRES_BUF
        r = crypto_sign_open(sm, &mlen, sm, smlen, pk, ver_buf);
        r |= crypto_sign_verify(sig, siglen, mi, i, pk, ver_buf);
#else 
        r = crypto_sign_open(sm, &mlen, sm, smlen, pk);
        r |= crypto_sign_verify(sig, siglen, mi, i, pk);
#endif
        if (r) {
            printf("ERROR: signature verification failed\n");
            return -1;
        }
        for (k = 0; k < i; k++) {
            if (sm[k] != mi[k]) {
                printf("ERROR: message recovery failed\n");
                return -1;
            }
        }
    }
    return 0;
}

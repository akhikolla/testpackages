#include "nb64.h"
#include "ctype.h"

namespace {
    const uint8_t hex[256] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        2, 3, 4, 5, 6, 7, 8, 9, 0, 0,
        0, 0, 0, 0, 0, 10, 11, 12, 13, 14,
        15, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 10, 11, 12,
        13, 14, 15, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0};
}

namespace DigitalNetNS {
    void hosttonb64(const uint64_t array[], int length, uint8_t byte[])
    {
        int k = 0;
        for (int i = 0; i < length; i++) {
            uint64_t tmp = array[i];
            for (int j = 1; j <= 8; j++) {
                uint8_t b = (tmp >> (64 - j * 8)) & 0xff;
                byte[k++] = b;
            }
        }
    }

    void nb64tohost(const uint8_t byte[], int length, uint64_t array[])
    {
        int k = 0;
        int len = length / 8;
        for (int i = 0; i < len; i++) {
            uint64_t tmp = 0;
            for (int j = 0; j < 8; j++) {
                tmp = tmp << 8;
                tmp = tmp | byte[k++];
            }
            array[i] = tmp;
        }
    }

    void hexchar64tohost(const uint8_t hexchar[], int length, uint64_t array[])
    {
        int k = 0;
        int len = length / 16;
        for (int i = 0; i < len; i++) {
            uint64_t tmp = 0;
            for (int j = 0; j < 8; j++) {
                uint8_t byte = 0;
                for (int z = 0; z < 2; z++) {
                    byte = byte << 4;
                    byte |= hex[hexchar[k++]];
                }
                tmp = tmp << 8;
                tmp = tmp | byte;
            }
            array[i] = tmp;
        }
    }

}

#pragma once
#ifndef NB64_H
#define NB64_H

#include <stdint.h>

namespace DigitalNetNS {
    void hosttonb64(uint64_t array[], int length, uint8_t byte[]);
    void nb64tohost(uint8_t byte[], int length, uint64_t array[]);
    void hexchar64tohost(const uint8_t hexchar[], int length, uint64_t array[]);
}
#endif // NB64_H

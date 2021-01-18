#ifndef _MD5_H
#define _MD5_H


#define uint8  unsigned char
#define uint32 unsigned long int

struct md5_context
{
    uint32 total[2];
    uint32 state[4];
    uint8 buffer[64];
};

void md5_starts( struct md5_context *ctx );
void md5_update( struct md5_context *ctx, uint8 *input, uint32 length );
void md5_finish( struct md5_context *ctx, uint8 digest[16] );

typedef struct md5_context MD5_CTX;


#endif /* md5.h */

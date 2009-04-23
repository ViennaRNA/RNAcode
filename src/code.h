#ifndef _CODE_H_
#define _CODE_H_

int* encodeSeq(char* seq);

int encodeAA(char aa);

char decodeAA(int encodedAA);

char* translateSeq(char* seq);

#endif

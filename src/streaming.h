#ifndef _STREAMING_H_
#define _STREAMING_H_

/**
 * Carries out the streaming step and writes the respective distribution functions from
 * collideField to streamField.
 */
void DoStreaming(float *collide_field, float *stream_field, int *flag_field, int xlength);
#endif

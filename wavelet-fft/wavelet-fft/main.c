//
//  main.c
//  wavelet-fft
//
//  Created by Helloworld Park. on 2017. 2. 25..
//  Copyright © 2017년 Helloworld Park. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <time.h>
#include "fourierTransform.h"

void compare(const unsigned long long q);

int main(int argc, const char * argv[]) {
    // insert code here...
    printf("Hello, World!\n");
    
    for (int i = 0; i < 100; i++) {
        for(unsigned long long q = 1; q <= 10; q++)
        {
            compare(1 << q);
        }
    }
    
    return 0;
}

void compare(const unsigned long long q)
{
    vecRef input = vecAlloc(q);
    for (unsigned long long i = 0; i < q; i++)
    {
        Real f = 0.0;
        f += 1.0 * sin((2.0 * M_PI * i) / q);
        f += 4.0 * cos((4.0 * M_PI * i) / q);
        f -= 2.0 * cos((16.0 * M_PI * i) / q);
        input->arr[i] = ComplexMake(f, 0);
    }
    vecRef output = NULL;
    clock_t tStart = clock();
    fft(input, &output);
    clock_t tEnd = clock();
    
    vecRef output2 = NULL;
    clock_t tStart2 = clock();
    dft(input, &output2);
    clock_t tEnd2 = clock();
    printf("N: %llu FFT: %lu DFT: %lu\n", q, tEnd - tStart, tEnd2 - tStart2);
    
    vecRelease(&input);
    vecRelease(&output);
    vecRelease(&output2);
}

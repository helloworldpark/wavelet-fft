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

typedef struct {
    clock_t fft;
    clock_t dft;
} elapsedtime;

elapsedtime compare(const unsigned long long q);
elapsedtime test01(const unsigned long long q);

int main(int argc, const char * argv[])
{
    for(unsigned long long q = 1; q <= 1; q++)
    {
        clock_t tFFT = 0.0;
        clock_t tDFT = 0.0;
        for (int i = 0; i < 1000; i++)
        {
            elapsedtime t = compare(1 << q);
            tFFT += t.fft;
            tDFT += t.dft;
        }
        printf("N: %llu FFT: %f DFT: %f\n", q, tFFT/1000.0, tDFT/1000.0);
    }
    
    test01(1<<8);
    
    return 0;
}

elapsedtime compare(const unsigned long long q)
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
    
    vecRelease(&input);
    vecRelease(&output);
    vecRelease(&output2);
    
    elapsedtime e;
    e.dft = tEnd2 - tStart2;
    e.fft = tEnd - tStart;
    return e;
}

elapsedtime test01(const unsigned long long q)
{
    vecRef input = vecAlloc(q);
    for (unsigned long long i = 0; i < q; i++)
    {
        Real f = 0.0;
        f += 0.3  * sin((2.0 * 1.0 * M_PI * i) / q);
        f += 4.2  * cos((2.0 * 5.0 * M_PI * i) / q);
        f -= 2.53 * cos((2.0 * 6.0 * M_PI * i) / q);
        input->arr[i] = ComplexMake(f, 0);
    }
    vecRef output1 = NULL;
    fft(input, &output1);
    vecRef output2 = NULL;
    ifft(output1, &output2);
    
    for(unsigned long long i = 0; i < output2->length; i++)
    {
        Real abs = cAbsolute(cSubtract(input->arr[i], output2->arr[i]));
        Complex b = input->arr[i];
        Complex c = output1->arr[i];
        printf("%llu %f %f %f %.10f\n", i, b.re, c.re, c.im, abs);
    }
    
    vecRelease(&input);
    vecRelease(&output1);
    
    elapsedtime e;
    return e;
}

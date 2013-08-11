DimensionalityReduction.jl is licensed under the MIT License:

> Copyright (c) 2012: John Myles White and other contributors.
>
> Permission is hereby granted, free of charge, to any person obtaining
> a copy of this software and associated documentation files (the
> "Software"), to deal in the Software without restriction, including
> without limitation the rights to use, copy, modify, merge, publish,
> distribute, sublicense, and/or sell copies of the Software, and to
> permit persons to whom the Software is furnished to do so, subject to
> the following conditions:
>
> The above copyright notice and this permission notice shall be
> included in all copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
> EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
> MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
> NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
> LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
> OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
> WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

The current implementation of NMF borrows several ideas from the
implementation in Luis Pedro Coelho's Milk package. That code is under
the MIT license as well.

The current implementation of ICA was translated from Vicente Zarzoso's
Matlab code. http://www.i3s.unice.fr/~zarzoso/robustica.html

REFERENCES:
- V. Zarzoso and P. Comon, <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/tnn10.pdf">"Robust independent component analysis by iterative maximization</a>
   <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/tnn10.pdf">of the kurtosis contrast with algebraic optimal step size"</a>,
  IEEE Transactions on Neural Networks, vol. 21, no. 2, pp. 248-261, Feb. 2010.
- V. Zarzoso and P. Comon, <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/ica07.pdf">"Comparative Speed Analysis of FastICA"</a>,
  in: Proceedings ICA-2007, 7th International Conference on Independent Component Analysis
      and Signal Separation, London, UK, September 9-12, 2007, pp. 293-300.
- V. Zarzoso, P. Comon and M. Kallel,  <a href = "http://www.i3s.unice.fr/~zarzoso/biblio/eusipco06.pdf">"How Fast is FastICA?"</a>,
  in: Proceedings EUSIPCO-2006, XIV European Signal Processing Conference,
      Florence, Italy, September 4-8, 2006.

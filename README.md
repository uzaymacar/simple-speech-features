# Simple 
Speech Features

## Introduction

`simple-speech-features` is a simple, straight-forward Python demonstration for extracting and manipulating 
acoustic and prosodic features from sound waves based on [`Praat`](http://www.fon.hum.uva.nl/praat/) 
and [`Parselmouth`](https://github.com/YannickJadoul/Parselmouth). All speech extraction utilities
and helper functions are included in `feature_extraction_utils.py`. The default values for parameters
and information in the docstrings are directly from `Praat`.

## Getting Started

1. Clone this repository: `git clone https://github.com/uzaymacar/speech-features.git`
2. Install Parselmouth: `pip install praat-parselmouth`
3. Start extracting speech features on your datasets as shown in `main.ipynb`!

## References
* Jadoul, Y., Thompson, B., & de Boer, B. (2018). Introducing Parselmouth: A Python interface to Praat. Journal of Phonetics, 71, 1-15. https://doi.org/10.1016/j.wocn.2018.07.001
* http://www.fon.hum.uva.nl/praat/manual/Query.html
* http://www.fon.hum.uva.nl/rob/NKI_TEVA/TEVA/HTML/Analysis.html
* Feinberg, D. R. (2019, October 8). Parselmouth Praat Scripts in Python. https://doi.org/10.17605/OSF.IO/6DWR3
* http://practicalcryptography.com/miscellaneous/machine-learning/guide-mel-frequency-cepstral-coefficients-mfccs/
* Sample sound file (`sample.wav`) acquired from http://www.voiptroubleshooter.com/open_speech/american.html

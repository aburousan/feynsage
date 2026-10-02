# feynsage against Kira 3.1

Same machine (ADRISHTA: AMD EPYC 9534, idle), same targets, same seeds, same masters.
Time = wall-clock of the whole run (feynsage includes about 2 s of Sage start-up), memory = peak resident set.
Every feynsage result was compared with Kira's, coefficient by coefficient.

## 1 thread

| family | r | s | targets | masters | feynsage ff | feynsage exact (trimmed) | Kira (Fermat) | Kira (FireFly) | identical |
|---|---|---|---|---|---|---|---|---|---|
| kite | 2 | 1 | 15 | 2 | 1.3 s / 347 MB | 1.3 s / 342 MB | 4.3 s / 75 MB | 2.2 s / 42 MB | yes |
| kite | 4 | 2 | 70 | 2 | 1.7 s / 351 MB | 1.4 s / 351 MB | 4.5 s / 82 MB | 2.4 s / 47 MB | yes |
| kite | 6 | 2 | 210 | 2 | 2.8 s / 372 MB | 1.6 s / 367 MB | 4.7 s / 82 MB | 2.9 s / 55 MB | yes |
| kite | 8 | 3 | 495 | 2 | 7.4 s / 429 MB | 2.3 s / 418 MB | 5.5 s / 98 MB | 4.4 s / 73 MB | yes |
| kite | 10 | 4 | 1001 | 2 | 17.5 s / 543 MB | 4.0 s / 526 MB | 8.0 s / 143 MB | 8.3 s / 145 MB | yes |
| vertex | 1 | 1 | 6 | 3 | 1.5 s / 351 MB | 1.5 s / 347 MB | 4.5 s / 156 MB | 2.4 s / 45 MB | yes |
| vertex | 2 | 2 | 21 | 3 | 2.4 s / 384 MB | 2.3 s / 366 MB | 5.1 s / 172 MB | 2.9 s / 60 MB | yes |
| vertex | 3 | 2 | 56 | 3 | 5.0 s / 422 MB | 3.5 s / 397 MB | 5.8 s / 189 MB | 3.5 s / 65 MB | yes |
| vertex | 4 | 3 | 126 | 3 | 15.1 s / 565 MB | 11.0 s / 540 MB | 9.5 s / 198 MB | 6.9 s / 151 MB | yes |
| vertex | 5 | 3 | 252 | 3 | 31.7 s / 714 MB | 23.1 s / 717 MB | 13.0 s / 222 MB | 9.8 s / 212 MB | yes |
| sunset | 2 | 1 | 12 | 3 | 1.3 s / 346 MB | 1.3 s / 344 MB | 4.3 s / 59 MB | 2.2 s / 36 MB | yes |
| sunset | 4 | 2 | 45 | 3 | 7.0 s / 355 MB | 415.7 s / 6227 MB | 4.4 s / 63 MB | 2.5 s / 44 MB | yes |
| sunset | 6 | 2 | 84 | 3 | 39.1 s / 384 MB | timeout | 4.6 s / 63 MB | 3.2 s / 68 MB | yes |
| sunset | 8 | 3 | 180 | 3 | 247.7 s / 466 MB | - | 5.7 s / 71 MB | 7.7 s / 141 MB | yes |
| sunset | 10 | 3 | 264 | 3 | 802.2 s / 587 MB | - | 7.7 s / 81 MB | 15.0 s / 233 MB | yes |
| doublebox | 0 | 1 | 2 | 1 | 1.7 s / 350 MB | 1.7 s / 347 MB | 5.2 s / 258 MB | 3.0 s / 62 MB | yes |
| doublebox | 1 | 1 | 14 | 8 | 6.3 s / 377 MB | timeout | 7.2 s / 328 MB | 4.9 s / 74 MB | yes |
| doublebox | 1 | 2 | 21 | 8 | 12.4 s / 420 MB | - | 13.5 s / 417 MB | 10.7 s / 136 MB | yes |
| doublebox | 2 | 2 | 84 | 8 | 58.1 s / 558 MB | - | 19.6 s / 452 MB | 16.1 s / 237 MB | yes |
| doublebox | 3 | 2 | 252 | 8 | 189.7 s / 799 MB | - | 34.8 s / 451 MB | 32.3 s / 384 MB | yes |

## 8 threads

| family | r | s | targets | masters | feynsage ff | feynsage exact (trimmed) | Kira (Fermat) | Kira (FireFly) | identical |
|---|---|---|---|---|---|---|---|---|---|
| kite | 2 | 1 | 15 | 2 | 1.4 s / 346 MB | - | 5.8 s / 105 MB | 2.9 s / 94 MB | yes |
| kite | 4 | 2 | 70 | 2 | 1.7 s / 354 MB | - | 5.8 s / 109 MB | 3.0 s / 99 MB | yes |
| kite | 6 | 2 | 210 | 2 | 2.6 s / 368 MB | - | 5.9 s / 121 MB | 3.2 s / 114 MB | yes |
| kite | 8 | 3 | 495 | 2 | 6.4 s / 420 MB | - | 6.3 s / 164 MB | 3.9 s / 163 MB | yes |
| kite | 10 | 4 | 1001 | 2 | 14.6 s / 534 MB | - | 7.5 s / 292 MB | 5.5 s / 290 MB | yes |
| vertex | 1 | 1 | 6 | 3 | 1.5 s / 349 MB | - | 5.9 s / 201 MB | 3.0 s / 96 MB | yes |
| vertex | 2 | 2 | 21 | 3 | 1.9 s / 367 MB | - | 6.1 s / 211 MB | 3.2 s / 114 MB | yes |
| vertex | 3 | 2 | 56 | 3 | 2.7 s / 396 MB | - | 6.3 s / 224 MB | 3.4 s / 146 MB | yes |
| vertex | 4 | 3 | 126 | 3 | 6.8 s / 532 MB | - | 7.5 s / 396 MB | 4.3 s / 426 MB | yes |
| vertex | 5 | 3 | 252 | 3 | 13.1 s / 700 MB | - | 8.6 s / 575 MB | 5.3 s / 573 MB | yes |
| sunset | 2 | 1 | 12 | 3 | 1.4 s / 346 MB | - | 5.8 s / 103 MB | 3.0 s / 99 MB | yes |
| sunset | 4 | 2 | 45 | 3 | 5.5 s / 357 MB | - | 5.8 s / 100 MB | 3.1 s / 103 MB | yes |
| sunset | 6 | 2 | 84 | 3 | 29.5 s / 381 MB | - | 5.9 s / 100 MB | 3.3 s / 109 MB | yes |
| sunset | 8 | 3 | 180 | 3 | 189.9 s / 464 MB | - | 6.4 s / 112 MB | 4.3 s / 191 MB | yes |
| sunset | 10 | 3 | 264 | 3 | 664.1 s / 574 MB | - | 7.3 s / 127 MB | 5.7 s / 286 MB | yes |
| doublebox | 0 | 1 | 2 | 1 | 1.8 s / 350 MB | - | 6.2 s / 299 MB | 3.2 s / 92 MB | yes |
| doublebox | 1 | 1 | 14 | 8 | 3.6 s / 365 MB | - | 6.7 s / 372 MB | 3.7 s / 142 MB | yes |
| doublebox | 1 | 2 | 21 | 8 | 5.6 s / 397 MB | - | 8.2 s / 449 MB | 5.1 s / 284 MB | yes |
| doublebox | 2 | 2 | 84 | 8 | 21.4 s / 582 MB | - | 9.9 s / 505 MB | 6.8 s / 473 MB | yes |
| doublebox | 3 | 2 | 252 | 8 | 70.4 s / 780 MB | - | 13.9 s / 745 MB | 11.4 s / 749 MB | yes |


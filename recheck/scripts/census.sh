#!/bin/sh
# class dimension and T_l factorisation mod l for all primes 13..199, with my vr2.py (mod-l class: U2 = 1/4, U3 = 1/9)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
for l in 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173 179 181 191 193 197 199; do
  python3 -I vr2.py $l 2>&1 | grep "class dim\|factor degrees\|repeated\|eigenspace"
done

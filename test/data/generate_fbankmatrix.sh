#!/bin/bash

if [[ -z $KALDI_ROOT ]]; then
  echo "Please set environment variable KALDI_ROOT to the root directory of Kaldi, which contains subdirectories: src, egs, ... of Kaldi"
  exit 1
fi

test_wav=en-us-hello.wav
fbank_text=/tmp/fbank_${test_wav}.txt
fbank_matrix=fbankmat_${test_wav}.txt
input_scp=/tmp/fbank-test.scp
config_file=/tmp/fbank-test.conf
output_ark=/tmp/fbank-test.ark

LF=$'\n'

echo "test-wav-1 ${test_wav}${LF}" > $input_scp

cat > $config_file << EOF
--num-mel-bins=40
--use-energy=false
--energy-floor=0.0
--raw-energy=true
--htk-compat=false
--use-log-fbank=true
--use-power=true
--sample-frequency=16000
--frame-length=25.0
--frame-shift=10.0
--preemphasis-coefficient=0.97
--dither=0.0  # Disable dither to make fbank matrix consistent
--window-type=hamming
--round-to-power-of-two=true
--snip-edges=true
--low-freq=20
--high-freq=8000
EOF

# Generate fbank feature using Kaldi 
$KALDI_ROOT/src/featbin/compute-fbank-feats --config=$config_file scp:$input_scp ark,scp:$output_ark,-
$KALDI_ROOT/src/featbin/copy-feats --binary=false ark:$output_ark ark,t:$fbank_text

# Extract fbank matrix from $fbank_text
if [[ -e $fbank_matrix ]]; then
  rm $fbank_matrix
fi
for t in $(cat $fbank_text); do
  [[ $t =~ ^[0-9\.]+(e[+-][0-9]+)?$ ]] && echo "${t}" >> $fbank_matrix
done

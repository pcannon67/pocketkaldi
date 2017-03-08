#!/bin/bash

test_wav_hello=$testdir/data/en-us-hello.wav
test_wav_cat=$testdir/data/en-us-cat.wav
config_file=$testdir/data/kaldi_fbank.conf
output_ark=fbank-2wav.ark
output_scp=fbank-2wav.scp
output_arktxt=fbank-2wav.ark.txt

input_scp=fbank-in2wav.scp
echo "test-wav-1 ${test_wav_hello}" > $input_scp
echo "test-wav-2 ${test_wav_cat}" >> $input_scp

# Extract the float-point values from kaldi ark file
function extract_float () {
    from_ark=$1
    to_text=$2

    ark_text=/tmp/ark_text.txt
    $kaldiroot/src/featbin/copy-feats ark:$from_ark ark,t:$ark_text
    if [[ -e $to_text ]]; then
      rm $to_text
    fi
    for t in $(cat $ark_text); do
      [[ $t =~ ^[0-9\.]+(e[+-][0-9]+)?$ ]] && echo "${t}" >> $to_text
    done
}

# Compute fbank feature using Kaldi
$kaldiroot/src/featbin/compute-fbank-feats --config=$config_file scp:$input_scp ark,scp:$output_ark,$output_scp
extract_float $output_ark ${output_arktxt}.1
mv $output_scp ${output_scp}.1

# Compute fbank feature using pocketkaldi
./compute_fbank scp:$input_scp ark,scp:$output_ark,$output_scp
extract_float $output_ark ${output_arktxt}.2
mv $output_scp ${output_scp}.2

diff ${output_scp}.1 ${output_scp}.2 > /dev/null 2>&1
[[ $? != 0 ]] && exit 232

# Check the ark file. Due to some float-point accuracy problem, we
# need to compare them in float-point comparasion way
paste ${output_arktxt}.1 ${output_arktxt}.2 | awk -F $'\t' '
function abs(v) {return v < 0 ? -v : v}
BEGIN {
    retcode = 0;
}
{
    if (abs($1 - $2) >= 0.0001) {
        retcode = 1;
        print $0
    }
}
END {
    exit retcode;
}'
[[ $? != 0 ]] && exit 233

# Success
exit 0

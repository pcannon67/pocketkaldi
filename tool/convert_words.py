import sys
import struct
import io

if len(sys.argv) != 3:
    print("Usage: python3 {}: <words.txt> <am-bin>".format(sys.argv[0]))
    print("Convert Kaldi words.txt to pocketkaldi format.")
    sys.exit(1)

from_file = sys.argv[1]
to_file = sys.argv[2]

words_idx = []
max_idx = 0
with open(from_file) as fd:
    for line in fd:
        fields = line.strip().split()
        if len(fields) != 2:
            print('2 columns expected, but {} found: {}'.format(len(fields), line))
        word, wid = fields
        wid = int(wid)
        words_idx.append((word, wid))
        if wid > max_idx:
            max_idx = wid
words = ['<empty>' for _ in range(max_idx + 1)]
for word, idx in words_idx:
    words[idx] = word

# Make the symbol buffer and symbol index for symboltable
fd_buffer = io.BytesIO()
symboltable_idx = []
symboltable_cntidx = 0
for word in words:
    word_bytes = (word.strip() + '\0').encode('utf-8')
    fd_buffer.write(word_bytes)
    symboltable_idx.append(symboltable_cntidx)
    symboltable_cntidx += len(word_bytes)
symboltable_buffer = fd_buffer.getbuffer()
assert(len(symboltable_buffer) == symboltable_cntidx)

# Dump the idx and buffer to file
with open(to_file, 'wb') as fd:
    fd.write(b"SYM0")
    fd.write(struct.pack("<i", 8 + len(symboltable_idx) * 4 + len(symboltable_buffer)))
    print('Section size = ', 8 + len(symboltable_idx) * 4 + len(symboltable_buffer))
    fd.write(struct.pack("<i", len(symboltable_idx)))
    fd.write(struct.pack("<i", len(symboltable_buffer)))
    for idx in symboltable_idx:
        fd.write(struct.pack("<i", idx))
    fd.write(symboltable_buffer)

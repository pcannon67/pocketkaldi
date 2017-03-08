import fst
import sys
import struct

if len(sys.argv) != 3:
    print('Usage: python3 {} <openfst-binfile> <output-binfile>'.format(sys.argv[0]))
    quit()

try:
    t = fst.read(sys.argv[1])
except:
    print('Unable to read openfst file: {}'.format(sys.argv[1]))

start_state = t.start
state_number = len(t)
finals = []
arcs = []
for state in t.states:
    if state.final:
        finals.append(float(state.final))
    else:
        finals.append(0)

    for arc in state.arcs:
        arcs.append((
            state.stateid,
            arc.nextstate,
            arc.ilabel,
            arc.olabel,
            float(arc.weight)))

arcs.sort()
state_arcidx = [-1] * len(t)
last_state = -1
for idx, arc in enumerate(arcs):
    state = arc[0]
    if state_arcidx[state] == -1:
        state_arcidx[state] = idx

print(state_number, state_arcidx[: 10])
assert(len(state_arcidx) == len(t) and len(finals) == len(t))


file_size = 0
with open(sys.argv[2], 'wb') as fd:
    # Magic number
    b = struct.pack("<i", 0x3323)
    file_size += len(b)
    fd.write(b)
    
    b = struct.pack("<i", state_number)
    file_size += len(b)
    fd.write(b)

    b = struct.pack("<i", start_state)
    file_size += len(b)
    fd.write(b)

    for final in finals:
        b = struct.pack("<f", final)
        file_size += len(b)
        fd.write(b)
    for arcidx in state_arcidx:
        b = struct.pack("<f", arcidx)
        file_size += len(b)
        fd.write(b)
    for arc in arcs:
        b = struct.pack("<iiif", arc[1], arc[2], arc[3], arc[4])
        file_size += len(b)
        fd.write(b)

    # Magic number
    b = struct.pack("<i", 0x3325)
    file_size += len(b)
    fd.write(b)

print(file_size)

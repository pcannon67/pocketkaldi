import fst
import sys
import struct


def print_usage():
    print('Usage: python3 {} <openfst-binfile> <output-binfile> [text|binary]'.format(sys.argv[0]))

if len(sys.argv) not in {3, 4}:
    print_usage()
    quit()
output_binary = True
if len(sys.argv) == 4:
    if sys.argv[3] == "text":
        output_binary = False
    elif sys.argv[3] == "binary":
        output_binary = True
    else:
        print_usage()
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
        finals.append(float("inf"))

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
assert(len(state_arcidx) == len(t) and len(finals) == len(t))

if output_binary:
    with open(sys.argv[2], 'wb') as fd:
        # Magic number
        fd.write(struct.pack("<i", 0x3323))
        
        fd.write(struct.pack("<i", state_number))
        fd.write(struct.pack("<i", len(arcs)))
        fd.write(struct.pack("<i", start_state))

        for final in finals:
            fd.write(struct.pack("<f", final))
        for arcidx in state_arcidx:
            fd.write(struct.pack("<i", arcidx))
        for arc in arcs:
            fd.write(struct.pack("<iiif", arc[1], arc[2], arc[3], arc[4]))
else:
    # With text format
    with open(sys.argv[2], 'w') as fd:
        fd.write("state_number = {}\n".format(state_number))
        fd.write("arc_number = {}\n".format(len(arcs)))
        fd.write("start_state = {}\n".format(start_state))

        fd.write("============ final =============\n")
        for idx, final in enumerate(finals):
            fd.write("{} -> {}\n".format(idx, final))
        fd.write("============ state_arcidx =============\n")
        for idx, arcidx in enumerate(state_arcidx):
            fd.write("{} -> {}\n".format(idx, arcidx))
        fd.write("============ arcs =============\n")
        for idx, arc in enumerate(arcs):
            fd.write("{} -> next_state({}), input_label({}), output_label({}), weight({})\n".format(idx, arc[1], arc[2], arc[3], arc[4]))

print("Success")

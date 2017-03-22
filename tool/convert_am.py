import sys
import re
import struct
import math

if len(sys.argv) != 3:
    print("Usage: python3 {}: <text-nnet2-am> <am-bin>".format(sys.argv[0]))
    print("Convert Kaldi nnet2 AM to pocketkaldi format.")
    print("    text-nnet2-am: The text format of nnet2 AM file could be obtained by kaldi/src/nnet2bin/nnet-am-copy.")
    sys.exit(1)

from_file = sys.argv[1]
to_file = sys.argv[2]


# Ids for different layers
LINEAR_LAYER = 0
RELU_LAYER = 1
NORMALIZE_LAYER = 2
SOFTMAX_LAYER = 3
ADD_LAYER = 4
MUL_LAYER = 5

class AM:
    def __init__(self):
        self.left_context = 0
        self.right_context = 0
        self.prior = None
        self.layers = []

    def layer_output_dim(self, layer, input_dim):
        layer_type = layer[0]
        if layer_type == LINEAR_LAYER:
            b = layer[2]
            return len(b)
        elif layer_type == RELU_LAYER:
            return input_dim
        elif layer_type == NORMALIZE_LAYER:
            return input_dim
        elif layer_type == SOFTMAX_LAYER:
            return input_dim
        elif layer_type == ADD_LAYER:
            return input_dim
        elif layer_type == MUL_LAYER:
            return input_dim

    def layer_input_dim(self, layer):
        layer_type = layer[0]
        if layer_type == LINEAR_LAYER:
            W = layer[1]
            return len(W[0])
        elif layer_type == MUL_LAYER:
            v = layer[1]
            return len(v)
        else:
            return None

    def verify(self):
        output_dim = None
        input_dim = None
        for idx, layer in enumerate(self.layers):
            expected_dim = self.layer_input_dim(layer)
            if output_dim != None and expected_dim != None and output_dim != expected_dim:
                raise Exception('input_dim == {} expected, but {} found in layer {}'.format(
                    expected_dim,
                    output_dim,
                    idx))
            input_dim = output_dim
            output_dim = self.layer_output_dim(layer, input_dim)

    def write_vector(self, fd, vec):
        fd.write(b"VEC0")
        fd.write(struct.pack("<i", len(vec) * 4 + 4))
        fd.write(struct.pack("<i", len(vec)))
        for v in vec:
            fd.write(struct.pack("<f", v))

    def write_matrix(self, fd, mat):
        fd.write(b"MAT0")
        fd.write(struct.pack("<i", 8))
        fd.write(struct.pack("<i", len(mat)))
        fd.write(struct.pack("<i", len(mat[0])))
        for col in mat:
            self.write_vector(fd, col)

    def write_layer(self, fd, layer):
        layer_type = layer[0]
        fd.write(b"LAY0")

        if layer_type == ADD_LAYER:
            scale = layer[1]
            fd.write(struct.pack("<i", 8))
            fd.write(struct.pack("<i", layer_type))
            fd.write(struct.pack("<f", scale))
        else:
            fd.write(struct.pack("<i", 4))
            fd.write(struct.pack("<i", layer_type))

        # Additional vectors and matrix
        if layer_type == LINEAR_LAYER:
            W = layer[1]
            b = layer[2]
            self.write_matrix(fd, W)
            self.write_vector(fd, b)
        elif layer_type == MUL_LAYER:
            v = layer[1]
            self.write_vector(fd, v)
        elif layer_type == ADD_LAYER:
            v = layer[2]
            self.write_vector(fd, v)

    def write(self, filename):
        with open(filename, 'wb') as fd:
            fd.write(b"AM~0")
            fd.write(struct.pack("<i", 8))
            fd.write(struct.pack("<i", self.left_context))
            fd.write(struct.pack("<i", self.right_context))
            fd.write(b"NNT0")
            fd.write(struct.pack("<i", 4))
            fd.write(struct.pack("<i", len(self.layers)))
            for layer in self.layers:
                self.write_layer(fd, layer)
            self.write_vector(fd, self.prior)

re_tag = re.compile(r'<(.*?)>(.*?)</(.*?)>', re.DOTALL)

def goto_token(token_name, text):
    re_token = re.compile(r'<{}>(.*)'.format(token_name), re.DOTALL)
    m = re_token.search(text)
    if m == None:
        raise Exception('unable to find token: {}'.format(token_name))
    return m.group(1)

def read_until_token(token_name, text):
    re_token = re.compile(r'(.*?)<{}>'.format(token_name), re.DOTALL)
    m = re_token.search(text)
    if m == None:
        raise Exception('unable to find token: {}'.format(token_name))
    return m.group(1)

def read_token(text):
    m = re.search(r'^\s*<(.*?)>(.*)', text, re.DOTALL)
    if m == None:
        raise Exception('read_token failed')
    return (m.group(1), m.group(2))

def read_int(text):
    m = re.search(r'^\s*(\d+)\s+(.*)', text, re.DOTALL)
    if m == None:
        raise Exception('read_int failed')
    return (int(m.group(1)), m.group(2))

def read_matrix(text, num_type = float):
    m = re.search(r'^\s*\[(.*?)\]\s*(.*)', text, re.DOTALL)
    if m == None:
        raise Exception('read_matrix failed')
    remained = m.group(2)
    text = m.group(1)
    lines = text.split('\n')
    matrix_cols = []
    row_num = 0
    for line in lines:
        if line.strip() == '': continue
        matrix_cols.append(list(map(num_type, line.strip().split())))
        if row_num == 0:
            row_num = len(matrix_cols[0])
        elif row_num != len(matrix_cols[-1]):
            raise Exception('Row number mismatch')
    return matrix_cols, remained

am = AM()

with open(from_file) as fd:
    model_text = fd.read()

# Nnet token
remained_text = goto_token('Nnet', model_text)
remained_text = read_until_token('/Nnet', remained_text)
remained_text = goto_token('NumComponents', remained_text)
num_components, remained_text = read_int(remained_text)
print('num_components = {}'.format(num_components))

# Components token
remained_text = goto_token('Components', remained_text)
remained_text = read_until_token('/Components', remained_text)

# Tokens in Components
while remained_text.strip() != '':
    token_tag, remained_text = read_token(remained_text)
    print(token_tag)
    end_tag = '/' + token_tag
    content_text = read_until_token(end_tag, remained_text)
    remained_text = goto_token(end_tag, remained_text)

    # Parse token_text
    if token_tag == 'SpliceComponent':
        content_text = goto_token('Context', content_text)
        context_mat, content_text = read_matrix(content_text, int)
        left_context = -context_mat[0][0]
        right_context = context_mat[0][-1]
        am.left_context = left_context
        am.right_context = right_context
        print('left_context = {}, right_context = {}'.format(left_context, right_context))
    elif token_tag == 'AffineComponentPreconditionedOnline':
        content_text = goto_token('LinearParams', content_text)
        W, content_text = read_matrix(content_text)
        print('W: {} * {}'.format(len(W), len(W[0])))
        content_text = goto_token('BiasParams', content_text)
        b, content_text = read_matrix(content_text)
        print('b: {} * {}'.format(len(b), len(b[0])))
        am.layers.append((LINEAR_LAYER, W, b[0]))
    elif token_tag == 'RectifiedLinearComponent':
        am.layers.append((RELU_LAYER, ))
    elif token_tag == 'NormalizeComponent':
        am.layers.append((NORMALIZE_LAYER, ))
    elif token_tag == 'FixedScaleComponent':
        content_text = goto_token('Scales', content_text)
        scales, content_text = read_matrix(content_text)
        print('Scales: {} * {}'.format(len(W), len(W[0])))
        am.layers.append((MUL_LAYER, scales[0]))
    elif token_tag == 'SoftmaxComponent':
        am.layers.append((SOFTMAX_LAYER, ))
    else:
        raise Exception('unexpected layer name: ' + token_tag)

# Prior
remained_text = goto_token('/Nnet', model_text)
prior, content_text = read_matrix(remained_text)
print('Prior: {} * {}'.format(len(prior), len(prior[0])))
am.prior = prior[0]

am.verify()
am.write(to_file)
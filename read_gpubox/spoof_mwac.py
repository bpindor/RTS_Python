from numpy import zeros, floor
n_ants = 128
npol = 2

from numpy import zeros

class map_t():
    def __init__(self):
        self.stn1 = -1
        self.stn2 = -1
        self.pol1 = -1
        self.pol2 = -1

def fill_mapping_matrix():

    corr_mapping = []
    for i in range(256):
        corr_mapping.append([None] * 256)

    single_pfb_mapping = []
    for i in range(64):
        single_pfb_mapping.append(int(floor(i/4)) + (i%4) * 16)

    pfb_output_to_input = zeros(256,dtype=int)
    for p in range(4):
        for inp in range(64):
            pfb_output_to_input[p*64 + inp] = (single_pfb_mapping[inp] + p*64)

    for inp1 in range(n_ants):
        for inp2 in range(n_ants):
            for pol1 in range(npol):
                for pol2 in range(npol):
                    index1 = inp1 * npol + pol1;
                    index2 = inp2 * npol + pol2;
                    corr_mapping[pfb_output_to_input[index1]][pfb_output_to_input[index2]] = map_t()
                    corr_mapping[pfb_output_to_input[index1]][pfb_output_to_input[index2]].stn1 = inp1
                    corr_mapping[pfb_output_to_input[index1]][pfb_output_to_input[index2]].stn2 = inp2
                    corr_mapping[pfb_output_to_input[index1]][pfb_output_to_input[index2]].pol1 = pol1
                    corr_mapping[pfb_output_to_input[index1]][pfb_output_to_input[index2]].pol2 = pol2

    return corr_mapping

def extract_full_matrix():

    matrix = {}
    for i in range(n_ants):
        for j in range(i+1):
            k = i * (i + 1) / 2 + j;
            for pol1 in range(npol):
                for pol2 in range(npol):
                    index = (k * npol + pol1) * npol + pol2;
                    matrix[(((0 * n_ants + i) * n_ants + j) * npol + pol1) * npol + pol2] = index
                    matrix[(((0 * n_ants + j) * n_ants + i) * npol + pol2) * npol + pol1] = index 

    return matrix

def load_rts_baselines(corr_mapping,matrix):

    rts2index = {}
    rts2dataIndex = {}

    for stn2_id in range(n_ants):
        for stn1_id in range(stn2_id):
            cc_ct = stn1_id + (stn2_id * (stn2_id-1)) / 2
            input_1 = stn1_id*2
            input_2 = stn2_id*2
            # consider p_prd = 0 only 
            the_mapping = corr_mapping[input_1][input_2]
            dataIndex = the_mapping.stn1 * n_ants*npol*npol + the_mapping.stn2 * npol*npol + the_mapping.pol1 * npol + the_mapping.pol2
            rts2index[cc_ct] = matrix[dataIndex]
            rts2dataIndex[cc_ct] = dataIndex

    return rts2index
                            
def load_auto_baselines(corr_mapping,matrix):

    rts2index = {}
    rts2dataIndex = {}

    for stn_id in range(n_ants):
        ac_ct = stn_id
        input_1 = stn_id*2
        input_2 = stn_id*2
        # consider p_prd = 0 only 
        the_mapping = corr_mapping[input_1][input_2]
        dataIndex = the_mapping.stn1 * n_ants*npol*npol + the_mapping.stn2 * npol*npol + the_mapping.pol1 * npol + the_mapping.pol2
        rts2index[ac_ct] = matrix[dataIndex]
        rts2dataIndex[ac_ct] = dataIndex

    return rts2index                            

import os

structure_dir = "structures/true_structures"
sequence_dir = "sequences"

def get_path_list(dir_name): 
    """
    """
    file_names = [os.path.splitext(name)[0] for name in os.listdir(dir_name)]
    return file_names

def change_name(dir_name, type_list, extension): 
    for name, rna_type in type_list: 
        new_name = (name.split('_'))
        new_name.insert(1, rna_type)
        new_name = "_".join(new_name)

        name = os.path.join(dir_name, name + extension)
        new_name = os.path.join(dir_name, new_name + extension)

        #print(name, new_name)
        
        os.rename(name, new_name)


    return

def change_type_list(type_list, type_map):
    type_list = [(t[0], type_map[t[1]]) for t in type_list]
    return type_list

def main(): 
    #file_names =  get_path_list(sequence_dir)

    #print(file_names)

    type_map = {0: "sncRNA", 1: "tRNA", 2: "16SrRNA", 3: "5SrRNA", 4: "other", 5: "otherrRNA"}

    type_list = [('104_bpRNARFAM8062', 4), 
                 ('104_bpRNASRP1', 0), 
                 ('105_bpRNACRW17910', 3),
                 ('106_URS00006A1431', 3),
                 ('106_URS0000A7D778', 3),
                 ('107_URS000063F1A1', 3),
                 ('108_bpRNAPDB403', 4), 
                 ('108_URS000067CB9E', 1), 
                 ('111_bpRNACRW18400', 3), 
                 ('111_bpRNACRW20144', 3), 
                 ('114_bpRNARFAM42796', 4),
                 ('117_bpRNARFAM1', 3),
                 ('124_bpRNACRW84955S', 3), 
                 ('134_URS000066550E', 3),
                 ('134_URS0000A066E6', 5),
                 ('151_bpRNAPDB499', 4),
                 ('152_URS000068D944', 0),
                 ('161_bpRNAPDB453WithPseudoknot', 4),
                 ('189_bpRNACRW10337WithPseudoknot', 2), 
                 ('189_URS0000733A5B', 5),
                 ('203_URS000074B037', 2),
                 ('223_bpRNARFAM38455', 4),
                 ('237_URS00006BA413', 0), 
                 ('259_bpRNACRW6632', 2), 
                 ('270_bpRNASRP79', 0), 
                 ('275_bpRNARNP123WithPseudoknot', 4), 
                 ('278_bpRNASRP45', 0), 
                 ('330_URS0000569A4A', 0), 
                 ('526_URS00007F6A16', 2),
                 ('653_URS0001002813', 5),
                 ('67_bpRNASPR1', 1),
                 ('67_URS00005D4795', 1),
                 ('69_URS00004A5357', 5), 
                 ('747_URS0001A7336C', 4), 
                 ('78_URS00003E17F2', 5), 
                 ('94_URS0000160754', 2)]
    
    type_list = change_type_list(type_list, type_map)
    
    change_name(structure_dir, type_list, ".dbn")
    change_name(sequence_dir, type_list, ".fasta")


if __name__ == '__main__': 
    main()
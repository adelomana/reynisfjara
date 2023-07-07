input_file = '/home/adrian/databases/ensembl/mouse_genome/Mus_musculus.GRCm39.dna.toplevel.fa'

chromosome_lines = 0
counting = False

with open(input_file, 'r') as f:
    for line in f:
        chromosome = line[1]

        if counting == True:
            chromosome_lines = chromosome_lines + 1
            print(line, chromosome_lines, chromosome_lines*60)

        if chromosome == '6':
            print(line)
            counting = True

        if 'GACATGCGGTGGAACAAG' in line:
            print(line)
            break

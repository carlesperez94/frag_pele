def main(pdb_path, lig_chain="L"):
    with open(pdb_path) as in_pdb:
        pdb = in_pdb.readlines()
    elements = []
    pdb_out = []
    dictionary_to_transcript = {}
    for line in pdb:
        if line[0:6] == "HETATM" and line[21:22] == lig_chain:
            element = line[76:78]
            elements.append(element)
            list_of_unique_elements = set(elements)
            counters = {i: elements.count(i) for i in list_of_unique_elements}
            counter = counters[line[76:78]]
            new_name = "{}{}".format(line[76:78].strip().upper(), counter)
            while True:
                new_name = new_name + " "
                if len(new_name) == 4:
                    break
            dictionary_to_transcript[line[12:16]] = new_name
            line = list(line)
            line[12:16] = new_name
            line = "".join(line)
            if len(new_name) > 4:
                raise ValueError("Length of the string {} is too long. Only 4 characters accepted.".format(new_name))
        pdb_out.append(line)

    return dictionary_to_transcript

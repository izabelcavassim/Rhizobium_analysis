# Creating a python code that takes the proteinortho ouput and transform it in fasta files
def read_proteinortho_groups(file_name):
    groups = []
    with open(file_name) as file_h:
        _ = file_h.__next__()

        print(_)

        for group_idx, l in enumerate(file_h):
            members_fields = l.strip().split("\t")[3:]
            members = []
            for member_field in members_fields:
                members.extend([(s, g) for s, g in (tuple(m.split("_")) for m in member_field.split(",") if m != "*")])

            groups.append((group_idx, members))

    return groups

bla = read_proteinortho_groups('/Users/PM/Downloads/myproject.proteinortho')
print(bla)

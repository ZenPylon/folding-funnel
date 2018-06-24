def get_amino_bonds(amino_name, start_index):
    """
    Construct a dictionary of bonds for the amino acid type specified.

    Args:
        amino_name: The name of the amino acid e.g. "GLY" 
        start_index: The number of atoms preceding this amino acid in the
                     polypeptide
    """
    amino_bonds = {}
    if amino_name == 'ILE':
        print('ILE')
    elif amino_name == 'GLY':
        print('GLY')
    else:
        print(f'ERROR - cannot create bonds for unsupported amino acid {amino_name}')


def get_amino_atom_count(amino_name):
    return _atom_counts[amino_name]

_atom_counts = {
    'GLY': 0,
    'ILE': 0
}
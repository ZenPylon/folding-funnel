import parmed as pmd
import parmed.charmm as charmm

# Create PSF from PDB file
filename = '1ubq'
ubiquitin = pmd.load_file(f'{filename}.pdb', structure=True)
ubiquitin.save(f'{filename}.psf', overwrite=True)

psf = charmm.CharmmPsfFile(f'{filename}.psf')

print(len(ubiquitin.atoms))
print(len(ubiquitin.angles))
print(len(ubiquitin.impropers))
print(len(psf.atoms))
print(len(psf.angles))
print(len(psf.impropers))


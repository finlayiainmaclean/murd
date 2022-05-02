import subprocess
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import PandasTools
import pandas as pd
import py3Dmol
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import prolif as plf
from prolif.plotting.network import LigNetwork
from src.utils import ncycle
def visualise_interactions(protein_file, lig_sdf_file):

    mol = Chem.MolFromPDBFile(protein_file, removeHs=True)
    prot = plf.Molecule(mol)
    cmd = f"obabel -isdf {lig_sdf_file} -opdb -O lig.pdb"
    subprocess.run(cmd, shell=True)
    mol = Chem.MolFromPDBFile("lig.pdb", removeHs=False)
    lig = plf.Molecule(mol)
    fp = plf.Fingerprint()
    fp.run_from_iterable([lig], prot)
    df = fp.to_dataframe(return_atoms=True)
    net = LigNetwork.from_ifp(df, lig,
                              # replace with `kind="frame", frame=0` for the other depiction
                              kind="frame", frame=0,
                              rotation=270)
    net.save("docked.html")
    interactions = [_1[1] + '//' + _1[2] for _1, _2 in fp.to_dataframe().to_dict().items()]

    return net.display(), interactions


def fix_protein(input_pdb_file, output_pdb_file):
    fixer = PDBFixer(filename=input_pdb_file)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb_file, 'w'))
    print("Finished!")


def prepare_receptor(prot_file):
    cmd = f"prepare_receptor4.py -r {prot_file}"
    subprocess.run(cmd, shell=True, check=True)
    print(prot_file)


def prepare_ligand(mol2_file):
    cmd = f"prepare_ligand4.py -l {mol2_file}"
    subprocess.run(cmd, shell=True, check=True)


def visualise_complex(prot_pdb_file=None, lig_sdf_files=None):
    v = py3Dmol.view()
    if prot_pdb_file is not None:
        v.addModel(open(prot_pdb_file).read())
        v.setStyle({'cartoon': {}, 'stick': {'radius': .1}})
        v.zoomTo({'model': 0})
    colors = ncycle(['greenCarbon', 'yellowCarbon', 'blueCarbon', 'redCarbon'])
    if lig_sdf_files is not None:
        for i, lig_sdf_file in enumerate(lig_sdf_files):
            v.addModelsAsFrames(open(lig_sdf_file).read())
            v.setStyle({'model': i + 1}, {'stick': {'colorscheme': next(colors)}})
        v.zoomTo({'model': 1})
    v.rotate(90)
    v.animate({'interval':1000})
    return v


class Dock:
    def __init__(self, prot_file, fix_pdb=False):
        self.output_sdf_file = None
        self.mol2_file = None
        self.smiles_file = None
        self.prot_file = prot_file
        if fix_pdb:
            fix_protein(prot_file, prot_file)
        prepare_receptor(prot_file)

    def dock(self, smiles, output_sdf_file='docked.sdf', com = None,
             size = (20,20,20), autobox_ligand = None, exhaustiveness=4,
             num_modes=9):

        self.smiles_file = 'lig.smiles'
        self.mol2_file = 'lig.mol2'
        self.output_sdf_file = output_sdf_file

        Path(self.smiles_file).write_text(smiles)

        cmd = f"""
        obabel \
        -ismi {self.smiles_file} \
        -omol2 \
        -O {self.mol2_file} \
        --gen3d --best --canonical \
        --conformers --weighted --nconf 50 --ff GAFF
        """
        subprocess.run(cmd, shell=True, check=True)

        prepare_ligand(self.mol2_file)

        sx, sy, sz = size
        cmd = f"""
        smina \
        -r {self.prot_file}qt -l {self.mol2_file.replace('.mol2', '.pdbqt')} \
        --num_modes {num_modes} -o {self.output_sdf_file} \
        --size_x {sx} --size_y {sy} --size_z {sz} \
        --exhaustiveness {exhaustiveness}
        """

        if com is not None:
            cx, cy, cz = com
            cmd += f"""--center_x {cx} --center_y {cy} --center_z {cz}"""
        elif autobox_ligand is not None:
            cmd += f"--autobox_ligand {autobox_ligand}"
        else:
            raise ValueError("Either use autobox_ligand or center_xyz")

        cmd = cmd.replace('\n','')
        print(cmd)
        subprocess.run(cmd, shell=True, check=True, capture_output=True)

        res = []
        for mol in Chem.SDMolSupplier(self.output_sdf_file):
            res.append((Chem.MolToSmiles(mol), mol.GetPropsAsDict()['minimizedAffinity']))

        output = pd.DataFrame(res)
        output.columns = ['smiles', 'affinity']
        PandasTools.AddMoleculeColumnToFrame(output, 'smiles', 'mol', includeFingerprints=False)

        return output

    def score(self, mol2_file, output_sdf_file='docked.sdf'):
        self.mol2_file = mol2_file
        self.output_sdf_file = output_sdf_file

        prepare_ligand(mol2_file)

        cmd = f"""
        smina \
        -r {self.prot_file} \
        -l {self.mol2_file.replace('.mol2', '.pdbqt')} \
        --score_only \
        --num_modes 1 \
        -o {self.output_sdf_file}
        """
        subprocess.run(cmd, shell=True, check=True, capture_output=True)

        with open(self.output_sdf_file, 'r') as f:
            lines = f.readlines()
        affinity = float(lines[-3])

        return affinity

    def visualise_complex(self):
        return visualise_complex(self.prot_file, [self.output_sdf_file])


    @classmethod
    def protein(cls, prot_file, **kwargs):
        return cls(prot_file, **kwargs)

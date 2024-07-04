import argparse
from ase import io
from ase.io import read
from ase.build import make_supercell
from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.bandstructure import HighSymmKpath

class SCFProcessor:
    def __init__(self, input_file, cif_file, mode='scf'):
        self.input_file = input_file
        self.cif_file = cif_file
        self.mode = mode
        self.filled_lines = []

    def ask_value(self, prompt):
        return input(prompt)
    
    def angstrom_to_bohr(self, angstrom):
       return angstrom * 1.8897259886

    def remove_cell_parameters(self, lines):
        cleaned_lines = []
        skip_lines = False
        for line in lines:
            if 'CELL_PARAMETERS angstrom' in line:
                skip_lines = True
                continue
            if skip_lines and line.strip() == '':
                skip_lines = False
                continue
            if not skip_lines:
                cleaned_lines.append(line)
        return cleaned_lines
    
    def remove_old_position(self, lines):
       cleaned_lines = []
       skip_lines = False
       for line in lines:
          if 'ATOMIC_POSITIONS angstrom' in line:
             cleaned_lines.append(line)  # Agrega la etiqueta pero activa el modo de salto
             skip_lines = True
             continue
          if skip_lines and line.strip() == '':  # Si se encuentra una línea en blanco, desactiva el modo de salto
             skip_lines = False
             continue
          if skip_lines and 'ATOMIC_POSITIONS' not in line:  # Salta líneas hasta encontrar otra etiqueta
              continue
          cleaned_lines.append(line)  # Agrega líneas no relacionadas con posiciones atómicas
       return cleaned_lines
    
    def get_kpoints_paths(self,cif_file):
      # Leer el archivo CIF y crear la estructura
      cif_parser = CifParser(cif_file)
      structure = cif_parser.get_structures()[0]

      # Usar HighSymmKpath para obtener los puntos de alta simetría y la ruta k
      kpath = HighSymmKpath(structure)
      kpoints = kpath.kpath

      # Extraer puntos de alta simetría y rutas
      kpoints_labels = kpoints['kpoints']
      kpaths = kpoints['path']

      # Crear una lista para almacenar las rutas k con etiquetas
      kpoints_paths = []

      for path in kpaths:
        for i in range(len(path) - 1):
            start_label = path[i]
            end_label = path[i + 1]
            start_point = kpoints_labels[start_label]
            end_point = kpoints_labels[end_label]
            kpoints_paths.append((start_label, start_point, end_label, end_point))

      return kpoints_paths


    def process_file(self):
        with open(self.input_file, 'r') as file:
            lines = file.readlines()

        inside_control_section = False
        inside_system_section = False
        inside_electrons_section = False
        inside_atomic_section = False
        inside_kpoints_section = False
        inside_atomP_section = False
        lines = self.remove_cell_parameters(lines)
        lines = self.remove_old_position(lines)
        for line in lines:
            if "&CONTROL" in line:
                inside_control_section = True
                self.filled_lines.append(line)
                self.fill_control_section(self.mode)
            elif inside_control_section and "/" in line:
                inside_control_section = False
                self.filled_lines.append(line)
            elif "&SYSTEM" in line:
                inside_system_section = True
                self.filled_lines.append(line)
            elif inside_system_section:
                if "ibrav" in line:
                    continue  # Omitir la línea original ibrav=0
                elif "/" in line:
                    self.fill_system_section()
                    inside_system_section = False
                    self.filled_lines.append(line)
                else:
                    self.filled_lines.append(line)
            elif "&ELECTRONS" in line:
                inside_electrons_section = True
                self.filled_lines.append(line)
                self.fill_electrons_section()
            elif inside_electrons_section and "/" in line:
                inside_electrons_section = False
                self.filled_lines.append(line)
            elif "ATOMIC_SPECIES" in line:
                inside_atomic_section = True
                self.filled_lines.append(line)
            elif inside_atomic_section:
                if "None" in line:
                    self.fill_atomic_species_section(line)
                elif line.strip() == '':
                    inside_atomic_section = False
                    self.filled_lines.append(line)
                else:
                    self.filled_lines.append(line)
            elif "K_POINTS gamma" in line:
                inside_kpoints_section = True
                self.fill_kpoints_section(self.mode)
            elif inside_kpoints_section and line.strip() == '':
                inside_kpoints_section = False
                self.filled_lines.append(line)
            elif "ATOMIC_POSITIONS angstrom" in line:
                inside_atomP_section = True
                self.filled_lines.append("ATOMIC_POSITIONS (crystal)\n")
                self.fill_ATOMIC_POSITIONS_section()
            elif inside_atomP_section and line.strip() == '':
                inside_atomP_section = False
                self.filled_lines.append(line)
            else:
                self.filled_lines.append(line)

        self.write_output()

    def fill_control_section(self, mode):
        if mode == "scf":
            self.filled_lines.append("   calculation      = 'scf'\n")
        if mode == "band":
            self.filled_lines.append("   calculation      = 'band'\n")
        pseudo_dir = self.ask_value("Please provide the value for pseudo_dir")
        self.filled_lines.append(f"   pseudo_dir       = '{pseudo_dir}'\n")
        self.filled_lines.append("   outdir           = 'temp'\n")
        prefix = self.ask_value("Please provide the value for prefix")
        self.filled_lines.append(f"   prefix           = '{prefix}'\n")
        self.filled_lines.append("   tefield          = .true.\n")

    def fill_system_section(self):
        ibrav = self.ask_value("Please provide the value for ibrav")
        self.filled_lines.append(f"   ibrav            = {ibrav}\n")
        structure = io.read(self.cif_file)
        cell = structure.get_cell()
        a = cell[0, 0]
        b = cell[1, 1]
        c = cell[2, 2]
        a_bohr = self.angstrom_to_bohr(a)
        b_bohr = self.angstrom_to_bohr(b)
        c_bohr = self.angstrom_to_bohr(c)
        if ibrav == "1" or ibrav == "2" or ibrav == "3":
            celldm1 = a_bohr
            self.filled_lines.append(f"   celldm(1)        = {celldm1}\n")
        if ibrav == "4" or ibrav == "6":
            celldm1 = a_bohr
            self.filled_lines.append(f"   celldm(1)        = {celldm1}\n")
            celldm3 = c_bohr / a_bohr
            self.filled_lines.append(f"   celldm(3)        = {celldm3}\n")
        if ibrav == "8":
            celldm1 = a_bohr
            self.filled_lines.append(f"   celldm(1)        = {celldm1}\n")
            celldm2 = b_bohr / a_bohr
            self.filled_lines.append(f"   celldm(2)        = {celldm2}\n")
            celldm3 = c_bohr / a_bohr
            self.filled_lines.append(f"   celldm(3)        = {celldm3}\n")
        self.filled_lines.append("   occupations      = 'smearing'\n")
        self.filled_lines.append("   smearing         = 'mv'\n")
        self.filled_lines.append("   degauss          = 0.02\n")
        ecutwfc = self.ask_value("Please provide the value for ecutwfc")
        self.filled_lines.append(f"   ecutwfc          = {ecutwfc}\n")
        edir = self.ask_value("Please provide the value for edir")
        self.filled_lines.append(f"   edir             = {edir}\n")
        self.filled_lines.append("   emaxpos          = 0\n")
        self.filled_lines.append("   eopreg           = 0.05\n")
        eamp = self.ask_value("Please provide the value for eamp")
        eamp2 = 0.0194469054 * float(eamp)
        self.filled_lines.append(f"   eamp             = {eamp2}\n")
        if self.mode == "band":
            nbnd = self.ask_value("Please provide the value for nbnd ")
            self.filled_lines.append(f"   nbnd             = {nbnd}\n")

    def fill_electrons_section(self):
        mixing_beta = self.ask_value("Please provide the value for mixing_beta")
        self.filled_lines.append(f"   mixing_beta      = {mixing_beta}\n")
        conv_thr = self.ask_value("Please provide the value for conv_thr")
        self.filled_lines.append(f"   conv_thr         = {conv_thr}\n")

    def fill_atomic_species_section(self, line):
        species, mass, _ = line.split()
        pseudoP = self.ask_value(f"Please provide the value PseudoP for {species} with mass {mass} ")
        self.filled_lines.append(f"{species} {mass} {pseudoP}\n")

    def fill_kpoints_section(self, mode):
        if mode == "scf":
          self.filled_lines.append("K_POINTS (automatic)\n")
          kpoint = self.ask_value("Please provide the value for k_point ")
          self.filled_lines.append(f"{kpoint} {kpoint} 1 0 0 0\n")
        if mode == "band":
          self.filled_lines.append("K_POINTS (crystal_b)\n")
          nPath = self.ask_value("Please provide the value for nPath")
          self.filled_lines.append(f"{nPath}\n")
          Lista = []
          points = []
          for i in range(int(nPath)):
              temp = self.ask_value(f"Please provide the value for nPath {i} ")
              temp2 = self.ask_value(f"Please provide the value for Point {i} ")
              Lista.append(f"{temp}")
              points.append(f"{temp2}")
          self.get_kpoints_paths(self.cif_file,Lista,points)
    
    def fill_ATOMIC_POSITIONS_section(self):
        # Obtener la celda y las coordenadas cartesianas
        structure = read(self.cif_file)
        cell = structure.get_cell()
        cartesian_positions = structure.get_positions()

        # Convertir coordenadas cartesianas a fraccionarias
        frac_positions = structure.get_scaled_positions(wrap=False)

        # Mostrar las coordenadas fraccionarias
        for i, atom in enumerate(structure):
            self.filled_lines.append(f"{atom.symbol} {frac_positions[i][0]} {frac_positions[i][1]} {frac_positions[i][2]}\n")
    
    def get_kpoints_paths(self, cif_file, Lista, points):
       # Leer el archivo CIF y crear la estructura
       cif_parser = CifParser(cif_file)
       structure = cif_parser.get_structures()[0]

       # Usar HighSymmKpath para obtener los puntos de alta simetría y la ruta k
       kpath = HighSymmKpath(structure)
       kpoints = kpath.kpath

       # Extraer puntos de alta simetría y rutas
       kpoints_labels = kpoints['kpoints']
       kpaths = kpoints['path']
       i = 0
       for path in Lista:
          self.filled_lines.append(f"{kpoints_labels[path][0]} {kpoints_labels[path][1]} {kpoints_labels[path][2]} {points[i]} ! {path}\n")
          i=i+1   

    
    def write_output(self):
        with open('output.in', 'w') as output_file:
            output_file.writelines(self.filled_lines)



def main():
    parser = argparse.ArgumentParser(description="SCF and Band file processor")
    parser.add_argument('-scf', metavar='input_file', type=str, help='Input file for SCF processing')
    parser.add_argument('-band', metavar='input_file', type=str, help='Input file for Band processing')
    parser.add_argument('-cif', metavar='input_file', type=str, help='Input file for .cif processing')
    
    args = parser.parse_args()

    if args.scf:
        processor = SCFProcessor(args.scf, args.cif, mode="scf")
        processor.process_file()
    elif args.band:
        processor = SCFProcessor(args.band, args.cif, mode="band")
        processor.process_file()
    else:
        print("Uso: python3 scf_file.py -scf archivo_de_entrada | -band archivo_de_entrada --mode [auto|user]")

if __name__ == "__main__":
    main()

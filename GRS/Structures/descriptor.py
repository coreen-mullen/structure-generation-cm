

from evaluators import Evaluators

def run_lammps(dgradflag):

            # simulation settings
            fname = file_prefix
            lmp.command("clear")
            lmp.command("info all out log")
            lmp.command('units  metal')
            lmp.command('atom_style  atomic')
            lmp.command("boundary    p p p")
            lmp.command("atom_modify    map hash")
            lmp.command('neighbor  2.3 bin')
            # boundary
            lmp.command('boundary  p p p')
            # read atoms
            lmp.command('read_data  %s.data' % fname )
            utypes = []
            for atom in atoms:
                if atom.symbol not in utypes:
                    utypes.append(atom.symbol)
            for ind,typ in enumerate(utypes):
                number = atomic_numbers[typ]
                mass = atomic_masses[number]
                lmp.command('mass   %d %f' % (ind+1,mass))

            lmp.command("pair_style     zero %f" % maxcut)
            lmp.command(f"pair_coeff     * *")

            if dgradflag:
                lmp.command(f"compute     pace all pace coupling_coefficients.yace 1 1")
            else:
                lmp.command(f"compute     pace all pace coupling_coefficients.yace 1 0")

            # run

            lmp.command(f"thermo         100")
            #lmp.command(f"run {nsteps}")
            lmp.command(f"run 0")

# declare compute pace variables

        dgradflag = 0
        run_lammps(dgradflag)
        lmp_pace = lmp.numpy.extract_compute("pace", LMP_STYLE_GLOBAL, LMP_TYPE_ARRAY)
        descriptor_grads = lmp_pace[ : len(atoms), : -1]
        return descriptor_grads
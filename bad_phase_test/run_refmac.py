import subprocess
import os 

def main(): 
    _args = []
    _args += ["HKLIN", "bad_phase_test/8ejo/8ejo_phases.mtz"]
    _args += ["XYZIN", "bad_phase_test/8ejo/8ejo_removed_na_water.pdb"]
    
    _args += ["HKLOUT", "./bad_phase_test/8ejo/8ejo_r1.mtz"]
    _args += ["XYZOUT", "./bad_phase_test/8ejo/8ejo_removed_na_water_r1.cif"]
    _args += ["XMLOUT", "./xmlout.xml"]
    labin = "FP=FP"
    labin += " SIGFP=SIGFP"
    labin += " FREE=FREE"
    labin += " PHIB=PHIC"
    labin += " FOM=FOM"
    _stdin = []
    _stdin.append("LABIN " + labin)
    _stdin.append(f"NCYCLES 1")
    _stdin.append("WEIGHT AUTO")
    _stdin.append("MAKE HYDR NO")
    _stdin.append("MAKE NEWLIGAND NOEXIT")
    _stdin.append("PHOUT")
    _stdin.append("PNAME modelcraft")
    _stdin.append("DNAME modelcraft")
    _stdin.append("END")

    with open("stdout.txt", "w") as out_stream:
        with open("stderr.txt", "w") as err_stream:
            process = subprocess.Popen(
                args=["/Applications/ccp4-8.0/bin/refmac5"] + _args,
                stdin=subprocess.PIPE if _stdin else None,
                stdout=out_stream,
                stderr=err_stream,
                encoding="utf8",
                env={**os.environ,},
                cwd=os.getcwd(),
                )
        if _stdin:
            for line in _stdin:
                process.stdin.write(line + "\n")
            process.stdin.close()
        process.wait()


if __name__ == "__main__":
    main()

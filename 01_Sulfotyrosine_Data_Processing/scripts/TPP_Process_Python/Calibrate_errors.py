
import sys
import statistics

def calibrate_errors(results_file):
    print("Calibrating", results_file)
    f = open(results_file, "r")
    rawfile_to_Da_shifts = {}


    counter =0
    isotope_peak = 1.003355
    for line in f:
        line = line.rstrip("\n")
        if counter != 0:
            cells=line.split("\t")
            rawfile = cells[0].split(".")[0]

            # da_shifts = []
            # if rawfile in rawfile_to_Da_shifts:
            #     da_shifts = rawfile_to_Da_shifts[rawfile]
            # da_shift = float(cells[5])
            # replaced code above with slightly more concise version
            da_shifts = rawfile_to_Da_shifts.get(rawfile, [])
            da_shift = float(cells[5])

            # TODO - calibration is not quite right for isotopes
            # if int(round(da_shift,0))==1:
            #     da_shift = da_shift - isotope_peak
            # if int(round(da_shift,0))==2:
            #     da_shift = da_shift -  (2 * isotope_peak)
            # if int(round(da_shift,0))==3:
            #     da_shift = da_shift -  (3 * isotope_peak)

            # if da_shift < 0.2:  #exclude those where isotope error has occured
            #     da_shifts.append(da_shift)
            # rawfile_to_Da_shifts[rawfile] = da_shifts

            # new code attempting to resolve.
            # Adjust for isotope peak errors
            # here we have picked quite an arbitrary 0.2 Da instead of rounding
            for i in range(1, 4):  # Adjust for 1, 2, and 3 Da shifts
                if abs(da_shift - i * isotope_peak) < 0.2:
                    da_shift -= i * isotope_peak
                    break

            # Filter out shifts greater than 0.2 (once again arbitrary)
            if abs(da_shift) < 0.2:
                da_shifts.append(da_shift)

            rawfile_to_Da_shifts[rawfile] = da_shifts

        counter += 1

    rawfile_to_median = {}
    for rawfile in rawfile_to_Da_shifts:
        da_shifts = rawfile_to_Da_shifts[rawfile]
        median_shift = statistics.median(da_shifts)
        rawfile_to_median[rawfile] = median_shift
        #print(rawfile,median_shift)
    f.close()

    f = open(results_file, "r")
    f_out = open(results_file[:-4]+"_calibrated.tsv", "w")
    counter=0

    for line in f:
        line = line.rstrip("\n")
        if counter != 0:
            cells=line.split("\t")
            rawfile = cells[0].split(".")[0]
            #print("error=",cells[5])
            da_error = float(cells[5])
            calibrated_da_error = da_error - rawfile_to_median[rawfile]
            # Apply additional filtering while writing
            if abs(calibrated_da_error) < 0.2:
                f_out.write(line + "\t" + str(calibrated_da_error) + "\n")
            # put the below line in the above conditional to get rid of values that slipped through our deisotoping loop
            # f_out.write(line + "\t" + str(calibrated_da_error) + "\n")
        else:
            f_out.write(line + "\t" + "calibrated_error" + "\n")
        counter += 1
    # close the files, good practice
    f.close()
    f_out.close()



if len(sys.argv)!= 2:
    print("Exit - expected usage\npython ",  sys.argv[0]," results.tsv\n")
else:
    calibrate_errors(sys.argv[1])
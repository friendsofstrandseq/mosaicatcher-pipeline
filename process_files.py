import os, sys, shutil

d = sys.argv[1]
for sample in os.listdir(d):
    if os.path.isdir(d + sample):
        for file in os.listdir(d + sample):
            # print(file)
            if file.endswith('.bam'):
                bam_file = file

                if os.path.islink(d + sample + '/all/' + bam_file):
                    print(sample)
                    try:
                        os.remove(d + sample + '/all/' + bam_file)
                        shutil.copy(d + sample + '/' + bam_file, d + sample + '/all/' + bam_file)
                    except:
                        sys.exit(sample, " error")
                elif os.path.isfile(d + sample + '/all/' + bam_file) is False:
                    shutil.copy(d + sample + '/' + bam_file, d + sample + '/all/' + bam_file)         

                if os.path.islink(d + sample + '/selected/' + bam_file):
                    print(sample)
                    try:
                        os.remove(d + sample + '/selected/' + bam_file)
                        shutil.copy(d + sample + '/' + bam_file, d + sample + '/selected/' + bam_file)

                    except:
                        sys.exit(sample, " error")
                elif os.path.isfile(d + sample + '/selected/' + bam_file) is False:
                    shutil.copy(d + sample + '/' + bam_file, d + sample + '/selected/' + bam_file)                        

            if file.endswith('.bai'):
                bai_file = file

                if os.path.islink(d + sample + '/all/' + bai_file):
                    print(sample)
                    try:
                        os.remove(d + sample + '/all/' + bai_file)
                        shutil.copy(d + sample + '/' + bai_file, d + sample + '/all/' + bai_file)
                    except:
                        sys.exit(sample, " error")
                elif os.path.isfile(d + sample + '/all/' + bai_file) is False:
                    shutil.copy(d + sample + '/' + bai_file, d + sample + '/all/' + bai_file)


                if os.path.islink(d + sample + '/selected/' + bai_file):
                    print(sample)
                    try:
                        os.remove(d + sample + '/selected/' + bai_file)
                        shutil.copy(d + sample + '/' + bai_file, d + sample + '/selected/' + bai_file)
                    except:
                        sys.exit(sample, " error")
                elif os.path.isfile(d + sample + '/selected/' + bai_file) is False:
                    shutil.copy(d + sample + '/' + bai_file, d + sample + '/selected/' + bai_file)

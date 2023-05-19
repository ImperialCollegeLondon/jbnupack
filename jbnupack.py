##Jonathan Bostock
#7 May 2023


#Import modules
import os
import os.path as op
import csv
import nupack

#This is the main loop, all other functions are contained within
def main():

    #Get the location of the directory to process
    loc = ""
    #loc = input("Enter the location of the files to be analysed")
    if loc == "":
        loc = os.getcwd()

    files = [op.join(loc,f) for f in os.listdir(loc) if f[-4:] == ".npi"]

    for f in files:

        f_out = f[0:-4] + " out.txt"
        #Assume we must run the file
        run = True

        #Check if an output file exists
        if op.isfile(f_out):
            #If so, check if it is older than input file
            #If it is older, it will have a smaller .getmtime value
            #Therefoer if the out file
            if op.getmtime(f) < op.getmtime(f_out):
                run = False

        #Only run if we have to
        if run:
            file_raw = open(f,"r",newline="")

            file_info = parse_file(file_raw)

            file_raw.close()

            file_results = get_results(file_info)

            if file_results != None:
                save_results(file_results, f_out)

##Here we have a bunch of simple functions which process data
#This stops re-use of code
def make_model(file_info):
    return nupack.Model(
        material=file_info["Conditions"]["material"],
        celsius=file_info["Conditions"]["celsius"],
        ensemble=file_info["Conditions"]["ensemble"],
        sodium=file_info["Conditions"]["sodium"],
        magnesium=file_info["Conditions"]["magnesium"])

#Here we can expand notation in the "(3.4" form to the "(((...." form
def expand(string):

    new_string = ""

    while len(string) > 0:

        char = string[0]
        string = string[1:]
        n_string = "0"

        while len(string) > 0:
            if string[0] in "0123456789":
                n_string += string[0]
                string = string[1:]
            else:
                break

        #God this is awful I'm sorry Turing
        #You didn't invent computers for me to do this
        n=max(int(n_string),1)
        char_string = "".join([char]*n)
        new_string = "".join([new_string,char_string])


    return(new_string)

#Makes a name for a sequence based on the first and last three bases
def make_name(sequence):

    if len(sequence) < 6:
        return sequence
    else:
        return ".".join([sequence[0:3],str(len(sequence)-6),sequence[-4:-1]])

#Parses a row and determines if it is empty, a comment, or normal
def row_type(row):

    if row == []:
        return "empty"
    if len(row)==1 and row[0]=="":
        return "empty"
    if row[0] == "":
        return "tabbed"
    if row[0][0] == "#":
        return "comment"
    else:
        return "normal"

#A safe version of float(string)
def try_numeric(string):

    try:
        return float(string)
    except:
        return string

#This is to handle the design stack
def make_target_strands(shape, name):
    #Now sort out strands, it's a bit of a pain
    strand_texts = shape.split("+")
    strands = []
    for i, st in enumerate(strand_texts):
        index = str(i+1) + " " + name
        strand_seq = "".join(["N"]*(len(st)))
        domain = nupack.Domain(strand_seq,
                               name="Domain " + index)
        strands.append(
            nupack.TargetStrand([domain],
                                name="Strand " + index))
    return strands

##These are the longer functions
#Workhorse file parser
def parse_file(file):

    strands = {}
    test_tubes = {}
    target_domains = {}
    target_strands = {}
    designs = {}
    design_tubes = {}
    conditions={"material":     "dna",
                "celsius":      25,
                "ensemble":     "stacking",
                "sodium":       0.1,
                "magnesium":    0}

    file_matrix = [line.strip().split("\t") for line in file]

    for i, row in enumerate(file_matrix):
        rt = row_type(row)
        rl = row[0].lower()

        if "conditions" in rl:

            conditions = get_conditions(file_matrix[i+1:])
            continue

        if "target domain" in rl:

            target_domains = add_target_domain(file_matrix[i+1:],
                                               target_domains)
            continue

        if "target strand" in rl:

            target_strands = add_target_strand(file_matrix[i+1:],
                                               target_strands, target_domains)
            continue

        if "strand" in rl:

            strands = add_strand(file_matrix[i+1:], strands)
            continue

        if "test tube" in rl:

            test_tubes = add_test_tube(file_matrix[i+1:], test_tubes, strands)
            continue

        if "design tube" in rl or "target tube" in rl:

            design_tubes = add_design_tube(file_matrix[i+1:],
                                           design_tubes,
                                           designs)
            continue

        if "design" in rl or "target" in rl:

            designs = add_design(file_matrix[i+1:], designs, target_strands)
            continue

    results = {"Conditions":    conditions,
               "Strands":       strands,
               "Test Tubes":    test_tubes,
               "Designs":       designs,
               "Design Tubes":  design_tubes}

    return results

##The next few functions are all called directly by file_parser

#Pulls out a set of conditions
def get_conditions(fm):

    conditions={"material":     "dna",
                "celsius":      25,
                "ensemble":     "stacking",
                "sodium":       0.1,
                "magnesium":    0}

    for row in fm:
        rt = row_type(row)
        if rt == "empty":
            break
        elif rt == "normal":
            conditions[row[-2]] = try_numeric(row[-1].lower())

    return conditions
#Adds a strand to a set of strands and returns the set
def add_strand(fm, strands):

    info = {"sequence": "",
            "name":     "Sequence " + str(len(strands)+1)}

    for row in fm:
        rt = row_type(row)
        if rt == "empty":
            break
        elif rt == "normal":
            info[row[-2].lower()] = row[-1]

    strands[info["name"]] = nupack.Strand(info["sequence"],name=info["name"])

    return strands

#Adds a test tube, can handle named strands and directly-inputted sequences
def add_test_tube(fm, test_tubes, strands):

    info = {"name":             "Tube " + str(len(test_tubes)+1),
            "max complex size": 1,
            "strands":          {}}

    for row in fm:
        rt = row_type(row)
        if rt == "empty":
            break
        elif rt=="normal":
            if row[-2].lower() in info.keys():
                info[row[-2].lower()] = try_numeric(row[-1])
            else:
                #This means we're dealing with a strand
                try:
                    #This looks for a strand in the strands dictionary
                    #Which has the same name as row[2]
                    #If found, it adds this to our info
                    info["strands"][strands[row[-1]]] = try_numeric(row[-2])
                except:
                    #Otherwise, create a new strand and add it to the system
                    new_strand = (
                        nupack.Strand(row[-1],
                                      name="".join([info["name"],
                                                    " Strand ",
                                                    str(len(info["strands"])+1)]))
                    )
                    info["strands"][new_strand] = try_numeric(row[-2])

    test_tubes[info["name"]] = nupack.Tube(
        strands=info["strands"],
        complexes=nupack.SetSpec(max_size=int(info["max complex size"])),
        name=info["name"])

    return test_tubes

def add_target_domain(fm, target_domains):

    info = {"name":             "Target Domain " + str(len(target_domains)+1),
            "sequence":         "",
            "length":           1}

    for row in fm:
        rt = row_type(row)
        rl = row[0].lower()
        if rt == "empty":
            break
        elif rt == "normal":
            info[row[0]] = try_numeric(row[1])

    info["sequence"] = expand(info["sequence"])

    if info["sequence"] == "":
        info["sequence"] = "".join(["N"]*int(info["length"]))

    target_domains[info["name"]] = nupack.Domain(info["sequence"],
                                                 name=info["name"])

    return target_domains

def add_target_strand(fm, target_strands, target_domains):
    info = {"name":             "Target Strand " + str(len(target_strands)+1),
            "length":           1,
            "sequence":         "",
            "domains":          []}

    for row in fm:
        rt = row_type(row)
        rl = row[0].lower()
        if rt == "empty":
            break
        elif rt == "normal":
            if rl in info.keys():
                info[rl] = try_numeric(row[1])
            if rl =="domain":
                try:
                    #Check if we already have this domain
                    #First check if it is a reverse complement section
                    if row[1][0] == "~":
                        info["domains"].append(~target_domains[row[1][1:]])
                    elif row[1][-1] in "*'":
                        info["domains"].append(~target_domains[row[1][:-1]])
                    else:
                        info["domains"].append(target_domains[row[1]])
                except:
                    #If not, make one locally
                    name = info["name"] + " Domain "
                    name = name + str(len(info["domains"])+1)
                    info["domains"].append(nupack.Domain(
                        row[1],
                        name=name))

    if info["sequence"] == "":
        info["sequence"] = "".join(["N"]*info["length"])

    if info["domains"] == []:
        info["domains"].append(nupack.Domain(info["sequence"],
                                             info["name"]+ " Domain 1"))

    target_strands[info["name"]] = nupack.TargetStrand(
        info["domains"],
        name=info["name"])

    return target_strands

#Adds a new design (aka target complex)
def add_design(fm, designs, target_strands):

    info = {"name":     "Design " + str(len(designs)+1),
            "shape":    "",
            "strands":  []}

    for row in fm:
        rt=row_type(row)
        rl = row[0].lower()
        if rt =="empty":
            break
        elif rt=="normal":
            if rl == "name":
                info["name"] = try_numeric(row[1])
            elif rl == "shape":
                #This means we're dealing with the shape row
                #Have to expand manually
                sh = expand(row[-1])
                info["shape"] = sh
            elif rl=="strand":
                try:
                    info["strands"].append(target_strands[row[1]])
                except:
                    info["strands"] += make_target_strands(row[1],
                                                           info["name"])

    if info["strands"] == []:

        info["strands"] = make_target_strands(info["shape"],
                                              info["name"])


    designs[info["name"]] = nupack.TargetComplex(
        info["strands"],
        info["shape"],
        name=info["name"])

    return designs


#Adds a new design-containing test tube
def add_design_tube(fm, design_tubes, designs):

    info = {"name":                 "Design Tube " + str(len(design_tubes)+1),
            "max off target size":  3,
            "designs":              {}}

    for row in fm:
        rt=row_type(row)
        if rt == "empty":
            break
        elif rt == "normal":
            if row[0].lower() in info.keys():
                info[row[-2]]=row[-1]
            else:
                #We are dealing with a design here
                info["designs"][designs[row[-1]]] = try_numeric(row[-2])


    design_tubes[info["name"]] = nupack.TargetTube(
        on_targets=info["designs"],
        off_targets=nupack.SetSpec(max_size=info["max off target size"]),
        name=info["name"])

    return design_tubes

#Runs the results of a "File parser" 
def get_results(file_info):

    model = make_model(file_info)

    if len(file_info["Test Tubes"]) > 0:
        results = nupack.tube_analysis(
            tubes = [tt[1] for tt in file_info["Test Tubes"].items()],
            model=model)

        return results

    elif len(file_info["Design Tubes"]) > 0:
        design_spec = nupack.tube_design(
            tubes = [dt[1] for dt in file_info["Design Tubes"].items()],
            model=model)

        results = design_spec.run(trials=1)
        #for some reason this likes to return results as a list
        #But we can handle it B-)
        return results


    return None

#Pretty simple: saves results to a file called "file out.txt"
def save_results(results, f_out):

    #If we have a list, iterate
    if type(results) == list:
        for i, result in enumerate(results):
            if i == 0:
                #First one must be named f_out to prevent re-running
                name = f_out
            else:
                #After ones get a number
                name = f_out[-4:] + " " + str(i+1) + ".txt"
            #Save each result
            result.save_text(name)
    else:
        #If no list, just save
        results.save_text(f_out)

if __name__ == "__main__":
    main()

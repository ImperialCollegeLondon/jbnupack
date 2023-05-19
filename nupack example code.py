##5-May-2023
#Testing Nupack

#Import Nupack
#This syntax allows us to use nupack functions
#Without having to specify that they are from nupack
from nupack import *

#First step is to define a "Model"
#Models are objects in the programming sense
#A model contains all of the non-sequence information

#Here we will define five parameters:
#Material is obvious
#Ensemble is to do with internals
#Celcius is also self-explanatory
#Sodium refers to /all/ of the monovalent ions
#Magnesium just refers to Mg2+

model_1 = Model(material="dna",
                ensemble="stacking",
                celsius=25,
                sodium=0.1,
                magnesium=0.0)

#Material defaults to "rna"
#Ensemble defaults to "stacking"
#celcius defaults to 37
#sodium defaults to 1.0
#magnesium defaults to 0.0

#Keeping in mind defaults we could specify an identical model like this:

model_2 = Model(material="dna",
                celsius=25,
                sodium=0.1)

#Before we can analyse anything, we need some sequences
#This also requres us to declare objects
#The name variable must be assigned for this to work

strand_A = Strand("CATTCTCTCAGTGTATTA", name="A")
strand_B = Strand("TAATACA", name="B")
strand_C = Strand("CTGAGAGAATG", name="C")

#We can call the ".nt()" method of a strand to get its length
#We can also call .name to get the full name
#Note that .nt() requires brackets as it is a function
#But .name does not, since it's just a variable

#This is python-specific, feel free to ignore
#But the point is to output our strands' lengths and names
out = "Strand {name} is {length} bases long"

for strand in [strand_A, strand_B, strand_C]:
    print(out.format(name=strand.name,
                     length=strand.nt()))

#To get an output, one must specify a "test tube ensemble"
#Firstly, set up a dictionary of strands and concentrations

strand_concs_1 = {strand_A: 1e-6,
                  strand_B: 5e-7,
                  strand_C: 1.5e-6}

#Test tubes are created with "Tube"
#This needs a dictionary of strands and concentrations, and a name
#To deal with complexes, we must use "SetSpec"
#This creates a list of possible complexes

tube_1 = Tube(strands=strand_concs_1,
              name="Test tube 1",
              complexes=SetSpec(max_size=3))
#Finally, we can run an analysis job
#This requires a model, and a list of tubes
tube_list = [tube_1]

analysis_1 = tube_analysis(tubes=tube_list,
                           model=model_1)

#We can "print" it to the console
print(analysis_1)

#But a better answer is probably to save it
analysis_1.save_text("Example Analysis 1.txt")

#Now let's try and make a very simple hairpin, with the structure
#((((((((((.....)))))(((((+))))))))))

#First we have to define "domains"
#Each domain has a sequence and name
#N is used as the fully degenerate nucleotide code
d_lb = Domain("NNNNN",    name = "left base")
d_st = Domain("NNNNN",    name = "stem")
d_lp = Domain("NNNNN",    name ="short loop")
d_rb = Domain("NNNNN",    name = "right base")

#But you can also fix domains
d_t_lp = Domain("TTTTT",  name="T-loop")

#Now package these up into strands
s_hairpin = TargetStrand([d_lb, d_st, d_lp, ~d_st, d_rb],
                       name = "hairpin strand")
s_base = TargetStrand([~d_lb, ~d_rb],
                      name = "base strand")

t_hairpin = TargetStrand([d_lb, d_st, d_t_lp, ~d_st, d_rb],
                         name = "T hairpin strand")

#Specify the shape string
shape = "(10.5)5(5+)10"
shape_long = "((((((((((.....)))))(((((+))))))))))"

#Now we define a "target complex"
complex_1 = TargetComplex([s_hairpin, s_base],
                          shape,
                          name="Hairpin")

complex_2 = TargetComplex([t_hairpin, s_base],
                          shape_long,
                          name="T-Hairpin")

#Finally, put the complexes in target tubes
#Make dictionaries, like one would with strands
complexes_1 = {complex_1: 1e-6}
complexes_2 = {complex_2: 5e-7}

d_tube_1 = TargetTube(on_targets = complexes_1,
                      off_targets = SetSpec(max_size=3),
                      name="Design tube 1")

d_tube_2 = TargetTube(on_targets = complexes_2,
                      off_targets = SetSpec(max_size=3),
                      name="Design tube 2")

#Now we want to actually get our results
#Put our TargetTubes into a list
d_tubes = [d_tube_1, d_tube_2]

#Now we build a design specification
#This is another object
design_spec = tube_design(tubes = d_tubes,
                          model=model_1)

#We must call a the ".run()" of our spec
#This gives us our actual design results
design_results = design_spec.run(trials=1)

design_results[0].save_text("Example Design 1.txt")

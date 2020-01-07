from os import listdir
from os.path import isfile, join

def repeat_char(times, char):
    res = ""
    count = 0
    while count < times:
        res += char
        count += 1
    return res


files = [f for f in listdir(".") if isfile(join(".", f))]
#print files
yaml_files = [s for s in files if ".yaml" in s]
png_files = [s for s in files if ".png" in s]

standalone = []
matched = []
for y in yaml_files:
    prefix = y.split(".")[0]
    if prefix+".png" in png_files:
        matched.append(prefix)
    else:
        standalone.append(prefix)
examples = ".. _yaml-examples:\n\n"
examples += "Ascent Actions Examples\n"
examples +="=======================\n"
examples +="\n"
for i in matched:
    f = open(i+".yaml", "r")
    comment = f.readline()[1:]
    f.close()
    examples += comment
    examples += repeat_char(len(comment), "-") + "\n"
    examples +="\n"
    examples += "YAML actions:\n\n"
    examples += ".. literalinclude:: examples/" + i + ".yaml\n"
    examples +="\n"
    examples += "Resulting image:\n\n"
    examples += ".. image:: examples/" + i + ".png\n"
    examples +="\n"

for i in standalone:
    f = open(i+".yaml", "r")
    comment = f.readline()[1:]
    f.close()
    examples += comment
    examples += repeat_char(len(comment), "-") + "\n"
    examples +="\n"
    examples += "YAML actions:\n\n"
    examples += ".. literalinclude:: examples/" + i + ".yaml\n"
    examples +="\n"


examples_file = open("../Examples.rst", "w")
examples_file.write(examples)
examples_file.close()

echo "spurious signal, Z_[spur] (%), chi2 / nDof" > C2".txt"
python BackgroundModelling.py -catName C2 -myy E1 -mjj E1 -yamlFile ConfigSpuriousC2.yaml
python BackgroundModelling.py -catName C2 -myy E1 -mjj E2 -yamlFile ConfigSpuriousC2.yaml
python BackgroundModelling.py -catName C2 -myy E2 -mjj E1 -yamlFile ConfigSpuriousC2.yaml
python BackgroundModelling.py -catName C2 -myy E2 -mjj E2 -yamlFile ConfigSpuriousC2.yaml

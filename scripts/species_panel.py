"""Curated Phylogenetically Paired Extreme-vs-Control Vertebrate Panel (Stage 1).

44 unique wild vertebrate species forming 23 independent Extreme-vs-Control sister pairs
across 4 orthogonal stressor axes designed to deconvolve Hypoxia, Cold, and Heat under
the Fitness-Beats-Truth (FBT) / Energy Resistance framework:
  1. High_Altitude_Hypoxia (8 pairs: chronic hypobaric hypoxia + cold)
  2. Diving_Subterranean_Hypoxia (6 pairs: intermittent hypoxia / ischemia-reperfusion, NO altitude cold)
  3. Arctic_Cold (5 pairs: extreme cold / proton-leak thermogenesis under normoxia)
  4. Desert_Heat (4 pairs: extreme thermal/arid stress under normoxia)
"""

from __future__ import annotations

# Each entry: (scientific_name, fallback_names, common_name, taxonomic_clade, stressor_axis, role, sister_pair_ids, habitat_notes)
SPECIES_PANEL = [
    # --- AXIS 1: HIGH-ALTITUDE HYPOXIA (+ COLD) [8 pairs] ---
    ("Anser indicus", [], "Bar-headed goose", "Anseriformes", "High_Altitude_Hypoxia", "Extreme", ["ALT_GOOSE"], "Trans-Himalayan high-altitude flight (~7,000-9,000 m); canonical COX3 W116R benchmark."),
    ("Anser anser", [], "Greylag goose", "Anseriformes", "Control", "Control", ["ALT_GOOSE"], "Lowland congeneric sister control for Anser indicus."),

    ("Pseudopodoces humilis", [], "Tibetan ground tit", "Passeriformes", "High_Altitude_Hypoxia", "Extreme", ["ALT_TIT"], "Endemic Tibetan Plateau passerine (3,300-5,500 m)."),
    ("Parus major", [], "Great tit", "Passeriformes", "Control", "Control", ["ALT_TIT"], "Lowland temperate parid sister control for Pseudopodoces humilis."),

    ("Bos grunniens", ["Bos mutus"], "Yak", "Artiodactyla_Bovini", "High_Altitude_Hypoxia", "Extreme", ["ALT_YAK"], "Tibetan Plateau bovid (3,500-5,500 m)."),
    ("Bos taurus", [], "Domestic cattle", "Artiodactyla_Bovini", "Control", "Control", ["ALT_YAK"], "Lowland congeneric sister control for Bos grunniens."),

    ("Pantholops hodgsonii", [], "Tibetan antelope (Chiru)", "Artiodactyla_Caprinae", "High_Altitude_Hypoxia", "Extreme", ["ALT_CHIRU"], "High-altitude Tibetan endurance runner (4,000-5,500 m)."),
    ("Ovis aries", [], "Domestic sheep", "Artiodactyla_Caprinae", "Control", "Control", ["ALT_CHIRU", "HEAT_ORYX"], "Temperate bovid sister control for Tibetan antelope and Arabian oryx."),

    ("Equus kiang", [], "Tibetan wild ass (Kiang)", "Perissodactyla", "High_Altitude_Hypoxia", "Extreme", ["ALT_KIANG"], "Tibetan Plateau equid (4,000-5,000 m)."),
    ("Equus caballus", [], "Domestic horse", "Perissodactyla", "Control", "Control", ["ALT_KIANG"], "Lowland congeneric sister control for Equus kiang."),

    ("Ochotona curzoniae", [], "Plateau pika", "Lagomorpha", "High_Altitude_Hypoxia", "Extreme", ["ALT_PIKA"], "Tibetan Plateau non-hibernating hypoxia specialist (3,200-5,000 m)."),
    ("Oryctolagus cuniculus", [], "European rabbit", "Lagomorpha", "Control", "Control", ["ALT_PIKA"], "Lowland lagomorph sister control for Ochotona curzoniae."),

    ("Panthera uncia", ["Uncia uncia"], "Snow leopard", "Carnivora_Felidae", "High_Altitude_Hypoxia", "Extreme", ["ALT_LEOPARD"], "High-altitude montane felid (3,000-5,500 m)."),
    ("Panthera tigris", [], "Tiger", "Carnivora_Felidae", "Control", "Control", ["ALT_LEOPARD"], "Lowland congeneric sister control for Panthera uncia."),

    ("Rhinopithecus bieti", [], "Black snub-nosed monkey", "Primates", "High_Altitude_Hypoxia", "Extreme", ["ALT_MONKEY"], "Highest-altitude non-human primate (3,800-4,700 m, Yunnan/Tibet)."),
    ("Macaca mulatta", [], "Rhesus macaque", "Primates", "Control", "Control", ["ALT_MONKEY"], "Lowland cercopithecid sister control for Rhinopithecus bieti."),

    # --- AXIS 2: DIVING & SUBTERRANEAN HYPOXIA / ISCHEMIA-REPERFUSION (NO ALTITUDE COLD) [6 pairs] ---
    ("Physeter macrocephalus", ["Physeter catodon"], "Sperm whale", "Cetacea_Physeteridae", "Diving_Subterranean_Hypoxia", "Extreme", ["DIV_SPERM"], "Extreme deep diver (>2,000 m, >90 min breath-hold ischemia-reperfusion)."),
    ("Hippopotamus amphibius", [], "Hippopotamus", "Artiodactyla_Whippomorpha", "Control", "Control", ["DIV_SPERM"], "Closest semi-aquatic/terrestrial sister taxon to Cetacea."),

    ("Ziphius cavirostris", [], "Cuvier's beaked whale", "Cetacea_Ziphiidae", "Diving_Subterranean_Hypoxia", "Extreme", ["DIV_BEAKED"], "Deepest-diving mammal on Earth (2,992 m, >130 min dives)."),
    ("Tursiops truncatus", [], "Bottlenose dolphin", "Cetacea_Delphinidae", "Control", "Control", ["DIV_BEAKED"], "Shallow coastal odontocete control for Ziphius cavirostris."),

    ("Leptonychotes weddellii", [], "Weddell seal", "Carnivora_Phocidae", "Diving_Subterranean_Hypoxia", "Extreme", ["DIV_WEDDELL"], "Antarctic deep/prolonged diver (>600 m, >80 min apnea)."),
    ("Canis lupus", [], "Gray wolf", "Carnivora_Caniformia", "Control", "Control", ["DIV_WEDDELL"], "Terrestrial caniform control for diving pinnipeds."),

    ("Mirounga angustirostris", ["Mirounga leonina"], "Northern elephant seal", "Carnivora_Phocidae", "Diving_Subterranean_Hypoxia", "Extreme", ["DIV_ELEPHANTSEAL"], "Continuous deep diver (>1,500 m, extreme apnea/reperfusion)."),
    ("Phoca vitulina", [], "Harbor seal", "Carnivora_Phocidae", "Control", "Control", ["DIV_ELEPHANTSEAL"], "Shallow coastal phocid sister control for Mirounga angustirostris."),

    ("Heterocephalus glaber", [], "Naked mole-rat", "Rodentia_Hystricognathi", "Diving_Subterranean_Hypoxia", "Extreme", ["SUB_NMR"], "Subterranean hypoxia/anoxia & extreme-longevity champion (survives 18 min 0% O2)."),
    ("Cavia porcellus", [], "Guinea pig", "Rodentia_Hystricognathi", "Control", "Control", ["SUB_NMR"], "Surface-dwelling hystricognath sister control for Heterocephalus glaber."),

    ("Nannospalax galili", ["Spalax galili", "Nannospalax ehrenbergi", "Spalax ehrenbergi"], "Blind mole-rat", "Rodentia_Muroidea", "Diving_Subterranean_Hypoxia", "Extreme", ["SUB_BMR"], "Subterranean hypoxia/hypercapnia specialist (survives 3% O2)."),
    ("Mus musculus", [], "House mouse", "Rodentia_Muroidea", "Control", "Control", ["SUB_BMR"], "Surface-dwelling muroid sister control for Nannospalax galili."),

    # --- AXIS 3: ARCTIC / EXTREME COLD (NORMOXIA) [5 pairs] ---
    ("Ursus maritimus", [], "Polar bear", "Carnivora_Ursidae", "Arctic_Cold", "Extreme", ["COLD_BEAR"], "High-Arctic sea-ice cold & high-lipid metabolism specialist."),
    ("Ursus americanus", [], "American black bear", "Carnivora_Ursidae", "Control", "Control", ["COLD_BEAR"], "Temperate congeneric sister control for Ursus maritimus."),

    ("Vulpes lagopus", ["Alopex lagopus"], "Arctic fox", "Carnivora_Canidae", "Arctic_Cold", "Extreme", ["COLD_FOX"], "Extreme Arctic tundra cold specialist (-50 C)."),
    ("Vulpes vulpes", [], "Red fox", "Carnivora_Canidae", "Control", "Control", ["COLD_FOX", "HEAT_FOX"], "Temperate congeneric sister control for BOTH Arctic fox and Sahara Fennec fox."),

    ("Rangifer tarandus", [], "Reindeer (Caribou)", "Artiodactyla_Cervidae", "Arctic_Cold", "Extreme", ["COLD_REINDEER"], "Circumpolar Arctic cervid."),
    ("Odocoileus virginianus", [], "White-tailed deer", "Artiodactyla_Cervidae", "Control", "Control", ["COLD_REINDEER"], "Temperate cervid sister control for Rangifer tarandus."),

    ("Mammuthus primigenius", [], "Woolly mammoth", "Proboscidea", "Arctic_Cold", "Extreme", ["COLD_MAMMOTH"], "Pleistocene Arctic mammoth (complete RefSeq mitogenome)."),
    ("Elephas maximus", [], "Asian elephant", "Proboscidea", "Control", "Control", ["COLD_MAMMOTH"], "Tropical/subtropical sister taxon for Mammuthus primigenius."),

    ("Ovibos moschatus", [], "Muskox", "Artiodactyla_Caprinae", "Arctic_Cold", "Extreme", ["COLD_MUSKOX"], "High-Arctic tundra caprine specialist."),
    ("Capra hircus", [], "Domestic goat", "Artiodactyla_Caprinae", "Control", "Control", ["COLD_MUSKOX"], "Temperate caprine sister control for Ovibos moschatus."),

    # --- AXIS 4: DESERT / EXTREME HEAT (NORMOXIA) [4 pairs] ---
    ("Vulpes zerda", ["Fennecus zerda"], "Fennec fox", "Carnivora_Canidae", "Desert_Heat", "Extreme", ["HEAT_FOX"], "Sahara hyper-arid thermal specialist; congener of temperate V. vulpes & Arctic V. lagopus."),

    ("Oryx leucoryx", [], "Arabian oryx", "Artiodactyla_Bovidae", "Desert_Heat", "Extreme", ["HEAT_ORYX"], "Arabian desert extreme heat & dehydration specialist."),

    ("Camelus dromedarius", [], "Dromedary camel", "Artiodactyla_Camelidae", "Desert_Heat", "Extreme", ["HEAT_CAMEL"], "Arid desert thermal/metabolic specialist."),
    ("Sus scrofa", [], "Domestic pig", "Artiodactyla_Suidae", "Control", "Control", ["HEAT_CAMEL"], "Temperate non-desert artiodactyl control for Camelus dromedarius."),

    ("Jaculus jaculus", [], "Lesser Egyptian jerboa", "Rodentia_Dipodidae", "Desert_Heat", "Extreme", ["HEAT_JERBOA"], "Sahara/Arabian desert thermal rodent specialist."),
    ("Rattus norvegicus", [], "Brown rat", "Rodentia_Muroidea", "Control", "Control", ["HEAT_JERBOA"], "Temperate rodent control for Jaculus jaculus."),
]

# Explicit pair definitions: (pair_id, stressor_axis, extreme_species, control_species, clade_group)
SISTER_PAIRS = [
    # High-Altitude Hypoxia (8)
    ("ALT_GOOSE", "High_Altitude_Hypoxia", "Anser indicus", "Anser anser", "Anseriformes"),
    ("ALT_TIT", "High_Altitude_Hypoxia", "Pseudopodoces humilis", "Parus major", "Passeriformes"),
    ("ALT_YAK", "High_Altitude_Hypoxia", "Bos grunniens", "Bos taurus", "Bovini"),
    ("ALT_CHIRU", "High_Altitude_Hypoxia", "Pantholops hodgsonii", "Ovis aries", "Caprinae_Altitude"),
    ("ALT_KIANG", "High_Altitude_Hypoxia", "Equus kiang", "Equus caballus", "Equidae"),
    ("ALT_PIKA", "High_Altitude_Hypoxia", "Ochotona curzoniae", "Oryctolagus cuniculus", "Lagomorpha"),
    ("ALT_LEOPARD", "High_Altitude_Hypoxia", "Panthera uncia", "Panthera tigris", "Felidae"),
    ("ALT_MONKEY", "High_Altitude_Hypoxia", "Rhinopithecus bieti", "Macaca mulatta", "Primates"),
    # Diving & Subterranean Hypoxia (6)
    ("DIV_SPERM", "Diving_Subterranean_Hypoxia", "Physeter macrocephalus", "Hippopotamus amphibius", "Physeteridae"),
    ("DIV_BEAKED", "Diving_Subterranean_Hypoxia", "Ziphius cavirostris", "Tursiops truncatus", "Ziphiidae"),
    ("DIV_WEDDELL", "Diving_Subterranean_Hypoxia", "Leptonychotes weddellii", "Canis lupus", "Phocidae_Weddell"),
    ("DIV_ELEPHANTSEAL", "Diving_Subterranean_Hypoxia", "Mirounga angustirostris", "Phoca vitulina", "Phocidae_Mirounga"),
    ("SUB_NMR", "Diving_Subterranean_Hypoxia", "Heterocephalus glaber", "Cavia porcellus", "Hystricognathi"),
    ("SUB_BMR", "Diving_Subterranean_Hypoxia", "Nannospalax galili", "Mus musculus", "Spalacidae"),
    # Arctic Cold (5)
    ("COLD_BEAR", "Arctic_Cold", "Ursus maritimus", "Ursus americanus", "Ursidae"),
    ("COLD_FOX", "Arctic_Cold", "Vulpes lagopus", "Vulpes vulpes", "Canidae_Arctic"),
    ("COLD_REINDEER", "Arctic_Cold", "Rangifer tarandus", "Odocoileus virginianus", "Cervidae"),
    ("COLD_MAMMOTH", "Arctic_Cold", "Mammuthus primigenius", "Elephas maximus", "Proboscidea"),
    ("COLD_MUSKOX", "Arctic_Cold", "Ovibos moschatus", "Capra hircus", "Caprinae_Arctic"),
    # Desert Heat (4)
    ("HEAT_FOX", "Desert_Heat", "Vulpes zerda", "Vulpes vulpes", "Canidae_Desert"),
    ("HEAT_ORYX", "Desert_Heat", "Oryx leucoryx", "Ovis aries", "Hippotraginae"),
    ("HEAT_CAMEL", "Desert_Heat", "Camelus dromedarius", "Sus scrofa", "Camelidae"),
    ("HEAT_JERBOA", "Desert_Heat", "Jaculus jaculus", "Rattus norvegicus", "Dipodidae"),
]

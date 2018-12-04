SeaDNA 12S rRNA reference library coverage
================
Rupert A. Collins
06 August 2018

##### Methods and description

This document describes the current sampling for the 12S rRNA fish reference library for the SeaDNA project. The document is a dynamic knitr document and can be updated quickly using the Makefile in SeaDNA/scripts. A list of species from the UK was generated from three sources: GBIF, FishBase, and the Water Framework Directive list of transitional species. This list was filtered to identify synonyms and duplicates, and annotated with FishBase taxonomic classification and FishBase common names. Next a sub-list of "common" species was generated. These were species that we believe are likely to be encountered in eDNA surveys of inshore and transitional waters of the UK, and comprise most of the species in Henderson (2015). Most of the remaining are either introduced species, rarely encountered migrants, oceanic pelagics, or deep sea organisms.

To calculate coverage we used the Bristol SeaDNA tissue catalogue, and also performed a search of the GenBank database. Because of inconsistencies in how researchers annotate their GenBank submissions and the differing internal coverage of primer pairs for particular gene fragments, we performed a search requesting all mitochondrial DNA. Then we pulled out the ~170 bp Miya fragment from all the mtDNA using a hidden markov model. This enabled us to have greater confidence that useful sequences had not been missed. For the resulting sequences we then tabulated all their metadata from GenBank in order to allow us to curate a custom reference library according to various criteria (e.g. must have reference specimen or locality data). <!-- NEED TO ALSO SEARCH GENBANK UNDER SYNONYMS, CHANGE/CHECK STREAKED GURNARD, INVESTIGATE REFSEQ SEQUENCES FOR REPETITION -->

##### Results

The total number of UK species is estimated to be around 538. Currently we have access to tissue samples or GenBank sequence data for 75% of the 171 common species (Table 2), and 57% of the 367 rare species (Table 3). For the common species, median number count tissues/sequences per species is 8, while 30 (18%) are represented by only one tissue or sequence. Table 1 lists the species of the highest priority to source.

###### Table 1. Highest priority species (common species with no tissue sample or GenBank data).

| Family          | Scientific Name               | Common Name            |
|:----------------|:------------------------------|:-----------------------|
| Congridae       | *Conger conger*               | European conger        |
| Atherinidae     | *Atherina presbyter*          | Sand smelt             |
| Belonidae       | *Belone belone*               | Garfish                |
| Esocidae        | *Esox lucius*                 | Northern pike          |
| Lotidae         | *Ciliata septentrionalis*     | Northern rockling      |
| Lotidae         | *Gaidropsarus mediterraneus*  | Shore rockling         |
| Lotidae         | *Gaidropsarus vulgaris*       | Three-bearded rockling |
| Gobiesocidae    | *Apletodon dentatus*          | Small-headed clingfish |
| Gobiesocidae    | *Diplecogaster bimaculata*    | Two-spotted clingfish  |
| Gobiesocidae    | *Lepadogaster candolii*       | Connemarra clingfish   |
| Gobiesocidae    | *Lepadogaster lepadogaster*   | Shore clingfish        |
| Osmeridae       | *Osmerus eperlanus*           | European smelt         |
| Ammodytidae     | *Ammodytes marinus*           | Lesser sand-eel        |
| Ammodytidae     | *Ammodytes tobianus*          | Small sandeel          |
| Ammodytidae     | *Gymnammodytes semisquamatus* | Smooth sandeel         |
| Ammodytidae     | *Hyperoplus lanceolatus*      | Great sandeel          |
| Blenniidae      | *Blennius ocellaris*          | Butterfly blenny       |
| Blenniidae      | *Parablennius gattorugine*    | Tompot blenny          |
| Callionymidae   | *Callionymus reticulatus*     | Reticulated dragonet   |
| Gobiidae        | *Gobiusculus flavescens*      | Two-spotted goby       |
| Gobiidae        | *Pomatoschistus lozanoi*      | Lozano's goby          |
| Gobiidae        | *Pomatoschistus norvegicus*   | Norway goby            |
| Gobiidae        | *Thorogobius ephippiatus*     | Leopard-spotted goby   |
| Labridae        | *Labrus bergylta*             | Ballan wrasse          |
| Labridae        | *Symphodus bailloni*          | Baillon's wrasse       |
| Mullidae        | *Mullus barbatus*             | Red mullet             |
| Stichaeidae     | *Chirolophis ascanii*         | Yarrell's blenny       |
| Scophthalmidae  | *Scophthalmus rhombus*        | Brill                  |
| Scophthalmidae  | *Zeugopterus punctatus*       | Topknot                |
| Scophthalmidae  | *Zeugopterus regius*          | Eckstrom's topknot     |
| Agonidae        | *Agonus cataphractus*         | Hooknose               |
| Liparidae       | *Liparis liparis*             | Striped seasnail       |
| Liparidae       | *Liparis montagui*            | Montagus seasnail      |
| Triglidae       | *Trigloporus lastoviza*       | Streaked gurnard       |
| Syngnathidae    | *Entelurus aequoreus*         | Snake pipefish         |
| Syngnathidae    | *Hippocampus guttulatus*      | Long-snouted seahorse  |
| Syngnathidae    | *Hippocampus hippocampus*     | Short snouted seahorse |
| Syngnathidae    | *Nerophis ophidion*           | Straightnose pipefish  |
| Petromyzontidae | *Lampetra planeri*            | European brook lamprey |
| Scyliorhinidae  | *Scyliorhinus stellaris*      | Nursehound             |
| Triakidae       | *Galeorhinus galeus*          | Tope shark             |
| Dasyatidae      | *Dasyatis pastinaca*          | Common stingray        |

###### Table 2. All common UK species with counts for tissue samples collected and sequence data obtained from GenBank (number of individuals).

| Family          | Scientific Name               | Common Name              |  Tissue Count|  GenBank Count|
|:----------------|:------------------------------|:-------------------------|-------------:|--------------:|
| Anguillidae     | *Anguilla anguilla*           | European eel             |             4|             59|
| Congridae       | *Conger conger*               | European conger          |              |               |
| Atherinidae     | *Atherina boyeri*             | Big-scale sand smelt     |             6|              3|
| Atherinidae     | *Atherina presbyter*          | Sand smelt               |              |               |
| Belonidae       | *Belone belone*               | Garfish                  |              |               |
| Clupeidae       | *Alosa alosa*                 | Allis shad               |              |              2|
| Clupeidae       | *Clupea harengus*             | Atlantic herring         |             6|             47|
| Clupeidae       | *Sardina pilchardus*          | European pilchard        |             1|              5|
| Clupeidae       | *Sprattus sprattus*           | European sprat           |             9|              2|
| Engraulidae     | *Engraulis encrasicolus*      | European anchovy         |              |              1|
| Cyprinidae      | *Abramis brama*               | Freshwater bream         |              |              3|
| Cyprinidae      | *Alburnus alburnus*           | Bleak                    |              |              2|
| Cyprinidae      | *Barbus barbus*               | Barbel                   |              |              2|
| Cyprinidae      | *Carassius auratus*           | Goldfish                 |              |             59|
| Cyprinidae      | *Carassius carassius*         | Crucian carp             |              |              5|
| Cyprinidae      | *Cyprinus carpio*             | Common carp              |              |             16|
| Cyprinidae      | *Gobio gobio*                 | Gudgeon                  |              |              2|
| Cyprinidae      | *Leuciscus idus*              | Ide                      |              |              1|
| Cyprinidae      | *Leuciscus leuciscus*         | Common dace              |              |              1|
| Cyprinidae      | *Phoxinus phoxinus*           | Eurasian minnow          |              |              7|
| Cyprinidae      | *Pseudorasbora parva*         | Stone moroko             |              |              5|
| Cyprinidae      | *Rutilus rutilus*             | Roach                    |              |              2|
| Cyprinidae      | *Scardinius erythrophthalmus* | Rudd                     |              |              2|
| Cyprinidae      | *Squalius cephalus*           | Chub                     |              |              2|
| Cyprinidae      | *Tinca tinca*                 | Tench                    |              |              4|
| Nemacheilidae   | *Barbatula barbatula*         | Stone loach              |              |              3|
| Esocidae        | *Esox lucius*                 | Northern pike            |              |               |
| Gadidae         | *Gadiculus argenteus*         | Silvery pout             |              |              2|
| Gadidae         | *Gadus morhua*                | Atlantic cod             |             1|            155|
| Gadidae         | *Melanogrammus aeglefinus*    | Haddock                  |              |              8|
| Gadidae         | *Merlangius merlangus*        | Whiting                  |            10|              6|
| Gadidae         | *Micromesistius poutassou*    | Blue whiting             |              |              7|
| Gadidae         | *Pollachius pollachius*       | Pollack                  |              |              5|
| Gadidae         | *Pollachius virens*           | Saithe                   |              |              9|
| Gadidae         | *Trisopterus esmarkii*        | Norway pout              |              |              3|
| Gadidae         | *Trisopterus luscus*          | Pouting                  |             3|              1|
| Gadidae         | *Trisopterus minutus*         | Poor cod                 |             6|              3|
| Lotidae         | *Ciliata mustela*             | Fivebeard rockling       |             7|              2|
| Lotidae         | *Ciliata septentrionalis*     | Northern rockling        |              |               |
| Lotidae         | *Enchelyopus cimbrius*        | Fourbeard rockling       |             1|              2|
| Lotidae         | *Gaidropsarus mediterraneus*  | Shore rockling           |              |               |
| Lotidae         | *Gaidropsarus vulgaris*       | Three-bearded rockling   |              |               |
| Lotidae         | *Molva molva*                 | Ling                     |              |              3|
| Merlucciidae    | *Merluccius merluccius*       | European hake            |             1|              3|
| Gasterosteidae  | *Gasterosteus aculeatus*      | Three-spined stickleback |              |              7|
| Gasterosteidae  | *Pungitius pungitius*         | Ninespine stickleback    |             1|              4|
| Gasterosteidae  | *Spinachia spinachia*         | Sea stickleback          |             5|              2|
| Gobiesocidae    | *Apletodon dentatus*          | Small-headed clingfish   |              |               |
| Gobiesocidae    | *Diplecogaster bimaculata*    | Two-spotted clingfish    |              |               |
| Gobiesocidae    | *Lepadogaster candolii*       | Connemarra clingfish     |              |               |
| Gobiesocidae    | *Lepadogaster lepadogaster*   | Shore clingfish          |              |               |
| Lophiidae       | *Lophius piscatorius*         | Angler                   |              |              3|
| Mugilidae       | *Chelon labrosus*             | Thicklip grey mullet     |              |              6|
| Mugilidae       | *Liza aurata*                 | Golden grey mullet       |              |              4|
| Mugilidae       | *Liza ramada*                 | Thinlip grey mullet      |             3|               |
| Osmeridae       | *Osmerus eperlanus*           | European smelt           |              |               |
| Ammodytidae     | *Ammodytes marinus*           | Lesser sand-eel          |              |               |
| Ammodytidae     | *Ammodytes tobianus*          | Small sandeel            |              |               |
| Ammodytidae     | *Gymnammodytes semisquamatus* | Smooth sandeel           |              |               |
| Ammodytidae     | *Hyperoplus immaculatus*      | Greater sand-eel         |             1|               |
| Ammodytidae     | *Hyperoplus lanceolatus*      | Great sandeel            |              |               |
| Blenniidae      | *Blennius ocellaris*          | Butterfly blenny         |              |               |
| Blenniidae      | *Coryphoblennius galerita*    | Montagu's blenny         |            17|               |
| Blenniidae      | *Lipophrys pholis*            | Shanny                   |            33|               |
| Blenniidae      | *Parablennius gattorugine*    | Tompot blenny            |              |               |
| Callionymidae   | *Callionymus lyra*            | Dragonet                 |            11|               |
| Callionymidae   | *Callionymus maculatus*       |                          |            10|               |
| Callionymidae   | *Callionymus reticulatus*     | Reticulated dragonet     |              |               |
| Caproidae       | *Capros aper*                 | Boarfish                 |             1|              2|
| Carangidae      | *Trachurus trachurus*         | Atlantic horse mackerel  |             8|              4|
| Cepolidae       | *Cepola macrophthalma*        | Red bandfish             |             6|               |
| Gobiidae        | *Aphia minuta*                | Transparent goby         |            10|              1|
| Gobiidae        | *Gobius cobitis*              | Giant goby               |              |              1|
| Gobiidae        | *Gobiusculus flavescens*      | Two-spotted goby         |              |               |
| Gobiidae        | *Gobius niger*                | Black goby               |             9|              3|
| Gobiidae        | *Gobius paganellus*           | Rock goby                |            11|              1|
| Gobiidae        | *Lesueurigobius friesii*      | Fries's goby             |              |              1|
| Gobiidae        | *Pomatoschistus lozanoi*      | Lozano's goby            |              |               |
| Gobiidae        | *Pomatoschistus microps*      | Common goby              |             1|               |
| Gobiidae        | *Pomatoschistus minutus*      | Sand goby                |              |              1|
| Gobiidae        | *Pomatoschistus norvegicus*   | Norway goby              |              |               |
| Gobiidae        | *Pomatoschistus pictus*       | Painted goby             |             5|               |
| Gobiidae        | *Thorogobius ephippiatus*     | Leopard-spotted goby     |              |               |
| Labridae        | *Centrolabrus exoletus*       | Rock cook                |              |              1|
| Labridae        | *Ctenolabrus rupestris*       | Goldsinny-wrasse         |              |              1|
| Labridae        | *Labrus bergylta*             | Ballan wrasse            |              |               |
| Labridae        | *Labrus mixtus*               | Cuckoo wrasse            |              |              1|
| Labridae        | *Symphodus bailloni*          | Baillon's wrasse         |              |               |
| Labridae        | *Symphodus melops*            | Corkwing wrasse          |             4|              1|
| Moronidae       | *Dicentrarchus labrax*        | European seabass         |             7|              4|
| Mullidae        | *Mullus barbatus*             | Red mullet               |              |               |
| Mullidae        | *Mullus surmuletus*           | Surmullet                |             1|               |
| Percidae        | *Gymnocephalus cernua*        | Ruffe                    |              |              2|
| Percidae        | *Perca fluviatilis*           | European perch           |              |             15|
| Percidae        | *Sander lucioperca*           | Pike-perch               |              |              7|
| Pholidae        | *Pholis gunnellus*            | Rock gunnel              |            10|              1|
| Scombridae      | *Scomber scombrus*            | Atlantic mackerel        |             2|              3|
| Sparidae        | *Pagellus bogaraveo*          | Blackspot seabream       |             1|              2|
| Sparidae        | *Sparus aurata*               | Gilthead seabream        |              |              4|
| Sparidae        | *Spondyliosoma cantharus*     | Black seabream           |             1|               |
| Stichaeidae     | *Chirolophis ascanii*         | Yarrell's blenny         |              |               |
| Stichaeidae     | *Lumpenus lampretaeformis*    | Snakeblenny              |              |              1|
| Trachinidae     | *Echiichthys vipera*          | Lesser weever            |             1|               |
| Trachinidae     | *Trachinus draco*             | Greater weever           |              |              2|
| Zoarcidae       | *Zoarces viviparus*           | Eelpout                  |              |              1|
| Bothidae        | *Arnoglossus laterna*         | Mediterranean scaldfish  |             5|               |
| Pleuronectidae  | *Glyptocephalus cynoglossus*  | Witch flounder           |             1|               |
| Pleuronectidae  | *Hippoglossus hippoglossus*   | Atlantic halibut         |             1|              6|
| Pleuronectidae  | *Limanda limanda*             | Common dab               |              |              1|
| Pleuronectidae  | *Microstomus kitt*            | Lemon sole               |             1|               |
| Pleuronectidae  | *Platichthys flesus*          | European flounder        |             5|               |
| Pleuronectidae  | *Pleuronectes platessa*       | European plaice          |            10|               |
| Scophthalmidae  | *Lepidorhombus whiffiagonis*  | Megrim                   |             1|               |
| Scophthalmidae  | *Phrynorhombus norvegicus*    | Norwegian topknot        |             9|               |
| Scophthalmidae  | *Scophthalmus maximus*        | Turbot                   |             1|              2|
| Scophthalmidae  | *Scophthalmus rhombus*        | Brill                    |              |               |
| Scophthalmidae  | *Zeugopterus punctatus*       | Topknot                  |              |               |
| Scophthalmidae  | *Zeugopterus regius*          | Eckstrom's topknot       |              |               |
| Soleidae        | *Buglossidium luteum*         | Solenette                |             9|               |
| Soleidae        | *Microchirus variegatus*      | Thickback sole           |             5|               |
| Soleidae        | *Pegusa lascaris*             | Sand sole                |             1|               |
| Soleidae        | *Solea solea*                 | Common sole              |             6|              1|
| Salmonidae      | *Oncorhynchus mykiss*         | Rainbow trout            |              |              2|
| Salmonidae      | *Salmo salar*                 | Atlantic salmon          |             1|              5|
| Salmonidae      | *Salmo trutta*                | Sea trout                |              |             15|
| Salmonidae      | *Thymallus thymallus*         | Grayling                 |              |              3|
| Agonidae        | *Agonus cataphractus*         | Hooknose                 |              |               |
| Cottidae        | *Cottus gobio*                | Bullhead                 |              |              1|
| Cottidae        | *Micrenophrys lilljeborgii*   | Norway bullhead          |             1|               |
| Cottidae        | *Myoxocephalus scorpius*      | Shorthorn sculpin        |              |              1|
| Cottidae        | *Taurulus bubalis*            | Longspined bullhead      |            12|              1|
| Cyclopteridae   | *Cyclopterus lumpus*          | Lumpfish                 |              |              2|
| Liparidae       | *Liparis liparis*             | Striped seasnail         |              |               |
| Liparidae       | *Liparis montagui*            | Montagus seasnail        |              |               |
| Triglidae       | *Chelidonichthys cuculus*     | Red gurnard              |             9|               |
| Triglidae       | *Chelidonichthys lucerna*     | Tub gurnard              |             2|               |
| Triglidae       | *Eutrigla gurnardus*          | Grey gurnard             |             7|               |
| Triglidae       | *Trigloporus lastoviza*       | Streaked gurnard         |              |               |
| Siluridae       | *Silurus glanis*              | Wels catfish             |              |              2|
| Syngnathidae    | *Entelurus aequoreus*         | Snake pipefish           |              |               |
| Syngnathidae    | *Hippocampus guttulatus*      | Long-snouted seahorse    |              |               |
| Syngnathidae    | *Hippocampus hippocampus*     | Short snouted seahorse   |              |               |
| Syngnathidae    | *Nerophis lumbriciformis*     | Worm pipefish            |             6|               |
| Syngnathidae    | *Nerophis ophidion*           | Straightnose pipefish    |              |               |
| Syngnathidae    | *Syngnathus acus*             | Greater pipefish         |             1|               |
| Syngnathidae    | *Syngnathus rostellatus*      | Nilsson's pipefish       |             1|               |
| Syngnathidae    | *Syngnathus typhle*           | Broadnosed pipefish      |              |              2|
| Balistidae      | *Balistes capriscus*          | Grey triggerfish         |              |              2|
| Molidae         | *Mola mola*                   | Ocean sunfish            |              |              4|
| Zeidae          | *Zeus faber*                  | John dory                |             5|              4|
| Petromyzontidae | *Lampetra fluviatilis*        | River lamprey            |              |              2|
| Petromyzontidae | *Lampetra planeri*            | European brook lamprey   |              |               |
| Petromyzontidae | *Petromyzon marinus*          | Sea lamprey              |              |              3|
| Carcharhinidae  | *Prionace glauca*             | Blue shark               |              |              5|
| Scyliorhinidae  | *Scyliorhinus canicula*       | Lesser spotted dogfish   |             2|              1|
| Scyliorhinidae  | *Scyliorhinus stellaris*      | Nursehound               |              |               |
| Triakidae       | *Galeorhinus galeus*          | Tope shark               |              |               |
| Triakidae       | *Mustelus asterias*           | Starry smooth-hound      |             2|               |
| Alopiidae       | *Alopias vulpinus*            | Thresher                 |              |              4|
| Cetorhinidae    | *Cetorhinus maximus*          | Basking shark            |              |             24|
| Lamnidae        | *Lamna nasus*                 | Porbeagle                |              |              3|
| Dasyatidae      | *Dasyatis pastinaca*          | Common stingray          |              |               |
| Rajidae         | *Amblyraja radiata*           | Starry ray               |              |              4|
| Rajidae         | *Dipturus batis*              | Blue skate               |              |              5|
| Rajidae         | *Leucoraja naevus*            | Cuckoo ray               |              |              3|
| Rajidae         | *Raja brachyura*              | Blonde ray               |              |              1|
| Rajidae         | *Raja clavata*                | Thornback ray            |             4|              3|
| Rajidae         | *Raja microocellata*          | Small-eyed ray           |              |              1|
| Rajidae         | *Raja montagui*               | Spotted ray              |              |              1|
| Rajidae         | *Raja undulata*               | Undulate ray             |              |              1|
| Squalidae       | *Squalus acanthias*           | Picked dogfish           |              |              3|

###### Table 3. All other UK species (reported in UK waters, but not listed in common species).

| Family             | Scientific Name                | Common Name                      |  Tissue Count|  GenBank Count|
|:-------------------|:-------------------------------|:---------------------------------|-------------:|--------------:|
| Acipenseridae      | *Acipenser ruthenus*           | Sterlet sturgeon                 |              |              4|
| Acipenseridae      | *Acipenser sturio*             | Sturgeon                         |              |              3|
| Anguillidae        | *Anguilla bengalensis*         | Indian mottled eel               |              |              4|
| Derichthyidae      | *Nessorhamphus ingolfianus*    | Duckbill oceanic eel             |              |              5|
| Muraenidae         | *Muraena helena*               | Mediterranean moray              |              |              3|
| Nemichthyidae      | *Avocettina infans*            | Avocet snipe eel                 |              |              5|
| Nemichthyidae      | *Nemichthys scolopaceus*       | Slender snipe eel                |              |              7|
| Serrivomeridae     | *Serrivomer beanii*            | Stout sawpalate                  |              |              3|
| Synaphobranchidae  | *Histiobranchus bathybius*     | Deep-water arrowtooth eel        |              |              1|
| Synaphobranchidae  | *Ilyophis blachei*             |                                  |              |               |
| Synaphobranchidae  | *Simenchelys parasitica*       | Snubnosed eel                    |              |              7|
| Synaphobranchidae  | *Synaphobranchus kaupii*       | Kaup's arrowtooth eel            |              |              7|
| Alepisauridae      | *Alepisaurus ferox*            | Long snouted lancetfish          |              |              5|
| Bathysauridae      | *Bathysaurus ferox*            | Deep-sea lizardfish              |              |              2|
| Evermannellidae    | *Evermannella balbo*           | Balbo sabretooth                 |              |               |
| Ipnopidae          | *Bathypterois dubius*          | Mediterranean spiderfish         |              |               |
| Notosudidae        | *Scopelosaurus lepidus*        | Blackfin waryfish                |              |              1|
| Paralepididae      | *Arctozenus risso*             | Spotted barracudina              |              |              1|
| Paralepididae      | *Magnisudis atlantica*         | Duckbill barracudina             |              |              2|
| Paralepididae      | *Paralepis coregonoides*       | Sharpchin barracudina            |              |               |
| Belonidae          | *Belone svetovidovi*           |                                  |              |               |
| Exocoetidae        | *Cheilopogon heterurus*        | Mediterranean flyingfish         |              |              1|
| Exocoetidae        | *Hirundichthys rondeletii*     | Black wing flyingfish            |              |              1|
| Scomberesocidae    | *Scomberesox saurus*           | Atlantic saury                   |              |               |
| Anoplogastridae    | *Anoplogaster cornuta*         | Common fangtooth                 |              |              3|
| Berycidae          | *Beryx decadactylus*           | Alfonsino                        |              |              4|
| Diretmidae         | *Diretmus argenteus*           | Silver spinyfin                  |              |              3|
| Trachichthyidae    | *Gephyroberyx darwinii*        | Darwin's slimehead               |              |              1|
| Trachichthyidae    | *Hoplostethus atlanticus*      | Orange roughy                    |              |               |
| Trachichthyidae    | *Hoplostethus mediterraneus*   | Mediterranean slimehead          |              |               |
| Clupeidae          | *Alosa fallax*                 | Twaite shad                      |              |               |
| Cobitidae          | *Cobitis taenia*               | Spined loach                     |              |              1|
| Cobitidae          | *Misgurnus fossilis*           | Weatherfish                      |              |               |
| Cyprinidae         | *Blicca bjoerkna*              | White bream                      |              |              4|
| Cyprinidae         | *Leucaspius delineatus*        | Belica                           |              |              2|
| Cyprinidae         | *Pimephales promelas*          | Fathead minnow                   |              |              7|
| Cyprinidae         | *Rhodeus amarus*               |                                  |              |              2|
| Cyprinidae         | *Rhodeus sericeus*             | Bitterling                       |              |              3|
| Poeciliidae        | *Poecilia reticulata*          | Guppy                            |              |               |
| Gadidae            | *Gadiculus thori*              |                                  |              |               |
| Gadidae            | *Raniceps raninus*             | Tadpole fish                     |              |              4|
| Lotidae            | *Brosme brosme*                | Tusk                             |              |              2|
| Lotidae            | *Gaidropsarus argentatus*      | Arctic rockling                  |              |              3|
| Lotidae            | *Gaidropsarus macrophthalmus*  | Bigeye rockling                  |              |               |
| Lotidae            | *Lota lota*                    | Burbot                           |              |              8|
| Lotidae            | *Molva dypterygia*             | Blue ling                        |              |               |
| Lotidae            | *Molva macrophthalma*          | Spanish ling                     |              |               |
| Macrouridae        | *Coelorinchus caelorhincus*    | Hollowsnout grenadier            |              |              1|
| Macrouridae        | *Coelorinchus caudani*         |                                  |              |               |
| Macrouridae        | *Coelorinchus labiatus*        | Spearsnouted grenadier           |              |              1|
| Macrouridae        | *Coelorinchus occa*            | Swordsnout grenadier             |              |              1|
| Macrouridae        | *Coryphaenoides armatus*       | Abyssal grenadier                |              |              1|
| Macrouridae        | *Coryphaenoides brevibarbis*   |                                  |              |              1|
| Macrouridae        | *Coryphaenoides guentheri*     | Gunther's grenadier              |              |              2|
| Macrouridae        | *Coryphaenoides mediterraneus* | Mediterranean grenadier          |              |               |
| Macrouridae        | *Coryphaenoides rupestris*     | Roundnose grenadier              |              |              4|
| Macrouridae        | *Macrourus berglax*            | Roughhead grenadier              |              |              2|
| Macrouridae        | *Malacocephalus laevis*        | Softhead grenadier               |              |              1|
| Macrouridae        | *Nezumia aequalis*             | Common Atlantic grenadier        |              |              1|
| Macrouridae        | *Nezumia sclerorhynchus*       | Roughtip grenadier               |              |               |
| Macrouridae        | *Trachyrincus murrayi*         | Roughnose grenadier              |              |              4|
| Macrouridae        | *Trachyrincus scabrus*         | Roughsnout grenadier             |              |              1|
| Melanonidae        | *Melanonus zugmayeri*          | Arrowtail                        |              |              3|
| Moridae            | *Antimora rostrata*            | Blue antimora                    |              |              3|
| Moridae            | *Guttigadus latifrons*         |                                  |              |               |
| Moridae            | *Halargyreus johnsonii*        | Slender codling                  |              |              2|
| Moridae            | *Lepidion eques*               | North Atlantic codling           |              |              1|
| Moridae            | *Mora moro*                    | Common mora                      |              |               |
| Phycidae           | *Phycis blennoides*            | Greater forkbeard                |              |              1|
| Gasterosteidae     | *Pungitius laevis*             | Smoothtail ninespine stickleback |              |              4|
| Lampridae          | *Lampris guttatus*             | Opah                             |              |              6|
| Regalecidae        | *Regalecus glesne*             | King of herrings                 |              |              1|
| Trachipteridae     | *Trachipterus arcticus*        | Dealfish                         |              |               |
| Ceratiidae         | *Ceratias holboelli*           | Kroyer's deep-sea angler fish    |              |               |
| Ceratiidae         | *Cryptopsaras couesii*         | Triplewart seadevil              |              |              3|
| Himantolophidae    | *Himantolophus groenlandicus*  | Atlantic footballfish            |              |              2|
| Linophrynidae      | *Photocorynus spiniceps*       |                                  |              |               |
| Lophiidae          | *Lophius budegassa*            | Blackbellied angler              |              |               |
| Melanocetidae      | *Melanocetus johnsonii*        | Humpback anglerfish              |              |              2|
| Oneirodidae        | *Chaenophryne draco*           | Smooth dreamer                   |              |              1|
| Oneirodidae        | *Oneirodes eschrichtii*        | Bulbous dreamer                  |              |               |
| Mugilidae          | *Mugil cephalus*               | Flathead grey mullet             |              |             24|
| Myctophidae        | *Benthosema glaciale*          | Glacier lantern fish             |              |              1|
| Myctophidae        | *Bolinichthys supralateralis*  | Stubby lanternfish               |              |               |
| Myctophidae        | *Ceratoscopelus maderensis*    | Madeira lantern fish             |              |              1|
| Myctophidae        | *Ceratoscopelus warmingii*     | Warming's lantern fish           |              |              3|
| Myctophidae        | *Diaphus dumerilii*            |                                  |              |              1|
| Myctophidae        | *Diaphus metopoclampus*        | Spothead lantern fish            |              |               |
| Myctophidae        | *Diaphus rafinesquii*          | White-spotted lantern fish       |              |              1|
| Myctophidae        | *Electrona risso*              | Electric lantern fish            |              |              2|
| Myctophidae        | *Lampadena speculigera*        | Mirror lanternfish               |              |               |
| Myctophidae        | *Lampanyctus crocodilus*       | Jewel lanternfish                |              |              2|
| Myctophidae        | *Lampanyctus intricarius*      | Diamondcheek lanternfish         |              |              2|
| Myctophidae        | *Lampanyctus macdonaldi*       | Rakery beaconlamp                |              |              1|
| Myctophidae        | *Lampanyctus pusillus*         | Pygmy lanternfish                |              |              1|
| Myctophidae        | *Myctophum nitidulum*          | Pearly lanternfish               |              |              2|
| Myctophidae        | *Myctophum punctatum*          | Spotted lanternfish              |              |              2|
| Myctophidae        | *Notolychnus valdiviae*        | Topside lampfish                 |              |              4|
| Myctophidae        | *Notoscopelus elongatus*       |                                  |              |              1|
| Myctophidae        | *Notoscopelus kroyeri*         | Lancet fish                      |              |               |
| Myctophidae        | *Protomyctophum arcticum*      | Arctic telescope                 |              |              1|
| Myctophidae        | *Symbolophorus veranyi*        | Large-scale lantern fish         |              |               |
| Neoscopelidae      | *Neoscopelus macrolepidotus*   | Large-scaled lantern fish        |              |              2|
| Halosauridae       | *Halosauropsis macrochir*      | Abyssal halosaur                 |              |              3|
| Halosauridae       | *Halosaurus johnsonianus*      | Halosaur                         |              |               |
| Notacanthidae      | *Notacanthus bonaparte*        | Shortfin spiny eel               |              |               |
| Notacanthidae      | *Notacanthus chemnitzii*       | Snubnosed spiny eel              |              |              7|
| Notacanthidae      | *Polyacanthonotus challengeri* | Longnose tapirfish               |              |              3|
| Notacanthidae      | *Polyacanthonotus rissoanus*   | Smallmouth spiny eel             |              |              3|
| Bythitidae         | *Cataetyx laticeps*            |                                  |              |               |
| Carapidae          | *Echiodon drummondii*          |                                  |              |              1|
| Ophidiidae         | *Brotula barbata*              | Bearded brotula                  |              |               |
| Ophidiidae         | *Lamprogrammus shcherbachevi*  | Scaleline cusk                   |              |               |
| Ophidiidae         | *Ophidion barbatum*            | Snake blenny                     |              |               |
| Ophidiidae         | *Spectrunculus grandis*        | Pudgy cuskeel                    |              |              1|
| Parabrotulidae     | *Leucobrotula adipata*         |                                  |              |               |
| Parabrotulidae     | *Parabrotula plagiophthalma*   | False cusk                       |              |               |
| Alepocephalidae    | *Alepocephalus agassizii*      | Agassiz' slickhead               |              |              3|
| Alepocephalidae    | *Alepocephalus bairdii*        | Baird's slickhead                |              |              2|
| Alepocephalidae    | *Alepocephalus rostratus*      | Risso's smooth-head              |              |               |
| Alepocephalidae    | *Bajacalifornia megalops*      | Bigeye smooth-head               |              |              3|
| Alepocephalidae    | *Narcetes stomias*             | Blackhead salmon                 |              |              2|
| Alepocephalidae    | *Rouleina attrita*             | Softskin smooth-head             |              |               |
| Alepocephalidae    | *Rouleina maderensis*          | Madeiran smooth-head             |              |               |
| Alepocephalidae    | *Xenodermichthys copei*        | Bluntsnout smooth-head           |              |              1|
| Argentinidae       | *Argentina silus*              | Greater argentine                |              |              1|
| Argentinidae       | *Argentina sphyraena*          | Argentine                        |              |               |
| Bathylagidae       | *Bathylagichthys greyae*       | Grey's deepsea smelt             |              |              1|
| Bathylagidae       | *Bathylagus euryops*           | Goiter blacksmelt                |              |               |
| Microstomatidae    | *Nansenia groenlandica*        | Greenland argentine              |              |               |
| Opisthoproctidae   | *Bathylychnops exilis*         | Javelin spookfish                |              |              3|
| Opisthoproctidae   | *Dolichopteroides binocularis* |                                  |              |               |
| Opisthoproctidae   | *Dolichopteryx rostrata*       |                                  |              |               |
| Platytroctidae     | *Barbantus curvifrons*         | Palebelly searsid                |              |               |
| Platytroctidae     | *Holtbyrnia anomala*           | Bighead searsid                  |              |              1|
| Platytroctidae     | *Holtbyrnia macrops*           | Bigeye searsid                   |              |               |
| Platytroctidae     | *Maulisia mauli*               | Maul's searsid                   |              |              2|
| Platytroctidae     | *Normichthys operosus*         | Multipore searsid                |              |              2|
| Platytroctidae     | *Sagamichthys schnakenbecki*   | Schnakenbeck's searsid           |              |               |
| Platytroctidae     | *Searsia koefoedi*             | Koefoed's searsid                |              |               |
| Ammodytidae        | *Ammodytes americanus*         | American sand lance              |              |              2|
| Ammodytidae        | *Gymnammodytes cicerelus*      | Mediterranean sand eel           |              |               |
| Anarhichadidae     | *Anarhichas denticulatus*      | Northern wolffish                |              |             23|
| Anarhichadidae     | *Anarhichas lupus*             | Atlantic wolffish                |              |             88|
| Anarhichadidae     | *Anarhichas minor*             | Spotted wolffish                 |              |             20|
| Bramidae           | *Brama brama*                  | Atlantic pomfret                 |              |               |
| Bramidae           | *Pterycombus brama*            | Atlantic fanfish                 |              |               |
| Bramidae           | *Taractes asper*               | Rough pomfret                    |              |              2|
| Bramidae           | *Taractichthys longipinnis*    | Big-scale pomfret                |              |               |
| Callanthiidae      | *Callanthias ruber*            | Parrot seaperch                  |              |               |
| Carangidae         | *Campogramma glaycos*          | Vadigo                           |              |               |
| Carangidae         | *Naucrates ductor*             | Pilotfish                        |              |              2|
| Carangidae         | *Seriola dumerili*             | Greater amberjack                |              |              6|
| Carangidae         | *Trachinotus ovatus*           | Pompano                          |              |              3|
| Centrarchidae      | *Ambloplites rupestris*        | Rock bass                        |              |              4|
| Centrarchidae      | *Lepomis gibbosus*             | Pumpkinseed                      |              |              6|
| Centrolophidae     | *Centrolophus niger*           | Rudderfish                       |              |               |
| Centrolophidae     | *Hyperoglyphe perciformis*     | Barrelfish                       |              |               |
| Centrolophidae     | *Schedophilus medusophagus*    | Cornish blackfish                |              |               |
| Chiasmodontidae    | *Chiasmodon niger*             | Black swallower                  |              |               |
| Chiasmodontidae    | *Pseudoscopelus altipinnis*    |                                  |              |               |
| Cichlidae          | *Cichlasoma bimaculatum*       | Black acara                      |              |               |
| Cichlidae          | *Oreochromis niloticus*        | Nile tilapia                     |              |               |
| Cichlidae          | *Tilapia zillii*               | Redbelly tilapia                 |              |               |
| Echeneidae         | *Remora remora*                | Shark sucker                     |              |              2|
| Epigonidae         | *Epigonus telescopus*          | Black cardinal fish              |              |               |
| Gempylidae         | *Lepidocybium flavobrunneum*   | Escolar                          |              |              4|
| Gempylidae         | *Nesiarchus nasutus*           | Black gemfish                    |              |              2|
| Gempylidae         | *Ruvettus pretiosus*           | Oilfish                          |              |              3|
| Gobiidae           | *Awaous commersoni*            |                                  |              |               |
| Gobiidae           | *Buenia jeffreysii*            | Jeffrey's goby                   |              |               |
| Gobiidae           | *Crystallogobius linearis*     | Crystal goby                     |              |              1|
| Gobiidae           | *Gobius couchi*                | Couch's goby                     |              |               |
| Gobiidae           | *Gobius cruentatus*            | Red-mouthed goby                 |              |              1|
| Gobiidae           | *Gobius gasteveni*             | Steven's goby                    |              |               |
| Gobiidae           | *Lebetus guilleti*             | Guillet's goby                   |              |               |
| Gobiidae           | *Lebetus scorpioides*          | Diminutive goby                  |              |               |
| Howellidae         | *Howella brodiei*              | Pelagic basslet                  |              |              2|
| Istiophoridae      | *Istiophorus albicans*         | Atlantic sailfish                |              |              2|
| Istiophoridae      | *Kajikia albida*               | Atlantic white marlin            |              |              4|
| Labridae           | *Acantholabrus palloni*        | Scale-rayed wrasse               |              |              1|
| Labridae           | *Coris julis*                  | Mediterranean rainbow wrasse     |              |              1|
| Lethrinidae        | *Lethrinus nebulosus*          | Spangled emperor                 |              |              3|
| Luvaridae          | *Luvarus imperialis*           | Luvar                            |              |              4|
| Nomeidae           | *Cubiceps gracilis*            | Driftfish                        |              |               |
| Nomeidae           | *Psenes maculatus*             | Silver driftfish                 |              |              2|
| Polyprionidae      | *Polyprion americanus*         | Wreckfish                        |              |               |
| Sciaenidae         | *Argyrosomus regius*           | Meagre                           |              |               |
| Sciaenidae         | *Umbrina cirrosa*              | Shi drum                         |              |               |
| Scombridae         | *Auxis rochei*                 | Bullet tuna                      |              |              6|
| Scombridae         | *Auxis thazard*                | Frigate tuna                     |              |              4|
| Scombridae         | *Euthynnus alletteratus*       | Little tunny                     |              |              3|
| Scombridae         | *Katsuwonus pelamis*           | Skipjack tuna                    |              |              4|
| Scombridae         | *Orcynopsis unicolor*          | Plain bonito                     |              |               |
| Scombridae         | *Sarda sarda*                  | Atlantic bonito                  |              |               |
| Scombridae         | *Thunnus alalunga*             | Albacore                         |              |             43|
| Scombridae         | *Thunnus albacares*            | Yellowfin tuna                   |              |              5|
| Scombridae         | *Thunnus thynnus*              | Atlantic bluefin tuna            |              |              7|
| Serranidae         | *Acanthistius sebastoides*     | Koester                          |              |               |
| Serranidae         | *Serranus cabrilla*            | Comber                           |              |               |
| Sparidae           | *Boops boops*                  | Bogue                            |              |               |
| Sparidae           | *Dentex dentex*                | Common dentex                    |              |              2|
| Sparidae           | *Dentex maroccanus*            | Morocco dentex                   |              |               |
| Sparidae           | *Pagellus acarne*              | Axillary seabream                |              |              2|
| Sparidae           | *Pagellus erythrinus*          | Common pandora                   |              |              2|
| Sparidae           | *Pagrus pagrus*                | Red porgy                        |              |              3|
| Sparidae           | *Sarpa salpa*                  | Salema                           |              |               |
| Stichaeidae        | *Leptoclinus maculatus*        | Daubed shanny                    |              |              5|
| Stromateidae       | *Peprilus triacanthus*         | Atlantic butterfish              |              |              3|
| Trichiuridae       | *Aphanopus carbo*              | Black scabbardfish               |              |              3|
| Trichiuridae       | *Lepidopus caudatus*           | Silver scabbardfish              |              |              1|
| Trichiuridae       | *Trichiurus lepturus*          | Largehead hairtail               |              |              2|
| Tripterygiidae     | *Tripterygion delaisi*         | Black-faced blenny               |              |              1|
| Xiphiidae          | *Xiphias gladius*              | Swordfish                        |              |              3|
| Zoarcidae          | *Lycenchelys sarsii*           | Sar's wolf eel                   |              |               |
| Zoarcidae          | *Lycodes esmarkii*             | Greater eelpout                  |              |               |
| Zoarcidae          | *Lycodes eudipleurostictus*    | Doubleline eelpout               |              |               |
| Zoarcidae          | *Lycodes pallidus*             | Pale eelpout                     |              |               |
| Zoarcidae          | *Lycodes squamiventer*         | Scalebelly eelpout               |              |               |
| Zoarcidae          | *Lycodes terraenovae*          |                                  |              |              1|
| Zoarcidae          | *Lycodes vahlii*               | Vahl's eelpout                   |              |              1|
| Zoarcidae          | *Lycodonus flagellicauda*      |                                  |              |               |
| Zoarcidae          | *Melanostigma atlanticum*      | Atlantic soft pout               |              |              1|
| Bothidae           | *Arnoglossus imperialis*       | Imperial scaldfish               |              |               |
| Bothidae           | *Arnoglossus thori*            | Thor's scaldfish                 |              |               |
| Pleuronectidae     | *Hippoglossoides platessoides* | American plaice                  |              |              3|
| Pleuronectidae     | *Reinhardtius hippoglossoides* | Greenland halibut                |              |             10|
| Scophthalmidae     | *Lepidorhombus boscii*         | Four-spot megrim                 |              |               |
| Soleidae           | *Bathysolea profundicola*      | Deepwater sole                   |              |               |
| Soleidae           | *Microchirus azevia*           |                                  |              |               |
| Soleidae           | *Microchirus theophila*        | Bastard sole                     |              |               |
| Salmonidae         | *Coregonus albula*             | Vendace                          |              |               |
| Salmonidae         | *Coregonus autumnalis*         | Arctic cisco                     |              |              2|
| Salmonidae         | *Coregonus clupeoides*         | Powan                            |              |               |
| Salmonidae         | *Coregonus lavaretus*          | European whitefish               |              |             81|
| Salmonidae         | *Coregonus pennantii*          | Gwyniad                          |              |               |
| Salmonidae         | *Coregonus pollan*             | Irish pollan                     |              |               |
| Salmonidae         | *Coregonus stigmaticus*        | Schelly                          |              |               |
| Salmonidae         | *Coregonus vandesius*          |                                  |              |               |
| Salmonidae         | *Oncorhynchus gorbuscha*       | Pink salmon                      |              |              5|
| Salmonidae         | *Salmo ferox*                  |                                  |              |               |
| Salmonidae         | *Salmo nigripinnis*            |                                  |              |               |
| Salmonidae         | *Salmo stomachicus*            |                                  |              |               |
| Salmonidae         | *Salvelinus alpinus*           | Arctic char                      |              |               |
| Salmonidae         | *Salvelinus fontinalis*        | Brook trout                      |              |             10|
| Salmonidae         | *Salvelinus gracillimus*       |                                  |              |               |
| Salmonidae         | *Salvelinus inframundus*       |                                  |              |               |
| Salmonidae         | *Salvelinus killinensis*       |                                  |              |               |
| Salmonidae         | *Salvelinus lonsdalii*         |                                  |              |               |
| Salmonidae         | *Salvelinus mallochi*          |                                  |              |               |
| Salmonidae         | *Salvelinus maxillaris*        |                                  |              |               |
| Salmonidae         | *Salvelinus namaycush*         | Lake trout                       |              |              8|
| Salmonidae         | *Salvelinus obtusus*           |                                  |              |               |
| Salmonidae         | *Salvelinus perisii*           |                                  |              |               |
| Salmonidae         | *Salvelinus struanensis*       |                                  |              |               |
| Salmonidae         | *Salvelinus willoughbii*       | Char                             |              |               |
| Salmonidae         | *Salvelinus youngeri*          | Golden charr                     |              |               |
| Cottidae           | *Artediellus atlanticus*       | Atlantic hookear sculpin         |              |               |
| Cottidae           | *Cottus perifretum*            | Chabot fluviatile                |              |              2|
| Cottidae           | *Icelus bicornis*              | Twohorn sculpin                  |              |               |
| Cottidae           | *Triglops murrayi*             | Moustache sculpin                |              |              1|
| Dactylopteridae    | *Dactylopterus volitans*       | Flying gurnard                   |              |              2|
| Liparidae          | *Careproctus longipinnis*      | Longfin snailfish                |              |               |
| Liparidae          | *Careproctus reinhardti*       | Sea tadpole                      |              |              1|
| Liparidae          | *Paraliparis bathybius*        | Black seasnail                   |              |               |
| Liparidae          | *Paraliparis challengeri*      |                                  |              |               |
| Liparidae          | *Paraliparis hystrix*          |                                  |              |               |
| Liparidae          | *Paraliparis kreffti*          |                                  |              |               |
| Peristediidae      | *Peristedion cataphractum*     | African armoured searobin        |              |               |
| Psychrolutidae     | *Cottunculus microps*          | Polar sculpin                    |              |               |
| Psychrolutidae     | *Cottunculus thomsonii*        | Pallid sculpin                   |              |              3|
| Scorpaenidae       | *Scorpaena porcus*             | Black scorpionfish               |              |               |
| Scorpaenidae       | *Scorpaena scrofa*             | Red scorpionfish                 |              |               |
| Sebastidae         | *Helicolenus dactylopterus*    | Blackbelly rosefish              |              |              1|
| Sebastidae         | *Sebastes norvegicus*          | Golden redfish                   |              |              1|
| Sebastidae         | *Sebastes viviparus*           | Norway redfish                   |              |               |
| Sebastidae         | *Trachyscorpia cristulata*     | Atlantic thornyhead              |              |               |
| Sebastidae         | *Trachyscorpia cristulata*     | Spiny scorpionfish               |              |               |
| Triglidae          | *Chelidonichthys obscurus*     | Longfin gurnard                  |              |               |
| Triglidae          | *Trigla lyra*                  | Piper gurnard                    |              |               |
| Ictaluridae        | *Ameiurus melas*               | Black bullhead                   |              |              2|
| Ictaluridae        | *Ameiurus nebulosus*           | Brown bullhead                   |              |              6|
| Ictaluridae        | *Ictalurus punctatus*          | Channel catfish                  |              |              7|
| Melamphaidae       | *Melamphaes suborbitalis*      | Shoulderspine bigscale           |              |               |
| Melamphaidae       | *Poromitra nigriceps*          |                                  |              |               |
| Melamphaidae       | *Scopelogadus beanii*          | Bean's bigscale                  |              |               |
| Melamphaidae       | *Scopelogadus mizolepis*       | Ragged bigscale                  |              |              2|
| Gonostomatidae     | *Bonapartia pedaliota*         | Longray fangjaw                  |              |              2|
| Gonostomatidae     | *Cyclothone braueri*           | Garrick                          |              |              3|
| Gonostomatidae     | *Cyclothone microdon*          | Veiled anglemouth                |              |              1|
| Gonostomatidae     | *Cyclothone obscura*           | Hidden bristlemouth              |              |              2|
| Gonostomatidae     | *Sigmops bathyphilus*          | Spark anglemouth                 |              |               |
| Sternoptychidae    | *Argyropelecus hemigymnus*     | Half-naked hatchetfish           |              |              1|
| Sternoptychidae    | *Argyropelecus olfersii*       |                                  |              |               |
| Sternoptychidae    | *Maurolicus muelleri*          | Silvery lightfish                |              |              1|
| Sternoptychidae    | *Sternoptyx diaphana*          | Diaphanous hatchet fish          |              |              4|
| Stomiidae          | *Borostomias antarcticus*      | Snaggletooth                     |              |              1|
| Stomiidae          | *Chauliodus sloani*            | Sloane's viperfish               |              |              5|
| Stomiidae          | *Flagellostomias boureei*      | Longbarb dragonfish              |              |               |
| Stomiidae          | *Malacosteus niger*            | Stoplight loosejaw               |              |              4|
| Stomiidae          | *Melanostomias bartonbeani*    | Scaleless black dragonfish       |              |               |
| Stomiidae          | *Rhadinesthes decimus*         | Slender snaggletooth             |              |               |
| Stomiidae          | *Stomias boa*                  | Boa dragonfish                   |              |              1|
| Stomiidae          | *Stomias boa*                  |                                  |              |              1|
| Stomiidae          | *Trigonolampa miriceps*        | Threelight dragonfish            |              |               |
| Centriscidae       | *Macroramphosus scolopax*      | Longspine snipefish              |              |              2|
| Molidae            | *Ranzania laevis*              | Slender sunfish                  |              |              4|
| Tetraodontidae     | *Lagocephalus lagocephalus*    | Oceanic puffer                   |              |              3|
| Tetraodontidae     | *Sphoeroides pachygaster*      | Blunthead puffer                 |              |              4|
| Grammicolepididae  | *Grammicolepis brachiusculus*  | Thorny tinselfish                |              |              1|
| Oreosomatidae      | *Neocyttus helgae*             | False boarfish                   |              |               |
| Pseudotriakidae    | *Pseudotriakis microdon*       | False catshark                   |              |              3|
| Scyliorhinidae     | *Apristurus aphyodes*          |                                  |              |               |
| Scyliorhinidae     | *Apristurus laurussonii*       | Iceland catshark                 |              |               |
| Scyliorhinidae     | *Apristurus profundorum*       | Deep-water catshark              |              |               |
| Scyliorhinidae     | *Galeus melastomus*            | Blackmouth catshark              |              |               |
| Scyliorhinidae     | *Galeus murinus*               | Mouse catshark                   |              |               |
| Sphyrnidae         | *Sphyrna zygaena*              | Smooth hammerhead                |              |              6|
| Triakidae          | *Mustelus mustelus*            | Smooth-hound                     |              |               |
| Chlamydoselachidae | *Chlamydoselachus anguineus*   | Frilled shark                    |              |              4|
| Hexanchidae        | *Heptranchias perlo*           | Sharpnose sevengill shark        |              |              4|
| Hexanchidae        | *Hexanchus griseus*            | Bluntnose sixgill shark          |              |             43|
| Lamnidae           | *Carcharodon carcharias*       | Great white shark                |              |            133|
| Lamnidae           | *Isurus oxyrinchus*            | Shortfin mako                    |              |              5|
| Dasyatidae         | *Pteroplatytrygon violacea*    | Pelagic stingray                 |              |              3|
| Myliobatidae       | *Mobula mobular*               | Devil fish                       |              |              4|
| Myliobatidae       | *Myliobatis aquila*            | Common eagle ray                 |              |              1|
| Arhynchobatidae    | *Bathyraja brachyurops*        | Broadnose skate                  |              |               |
| Arhynchobatidae    | *Bathyraja richardsoni*        | Richardson's ray                 |              |              2|
| Arhynchobatidae    | *Bathyraja spinicauda*         | Spinytail skate                  |              |              1|
| Rajidae            | *Amblyraja hyperborea*         | Arctic skate                     |              |              1|
| Rajidae            | *Amblyraja jenseni*            | Shorttail skate                  |              |               |
| Rajidae            | *Dipturus nidarosiensis*       | Norwegian skate                  |              |              3|
| Rajidae            | *Dipturus oxyrinchus*          | Longnosed skate                  |              |              4|
| Rajidae            | *Leucoraja circularis*         | Sandy ray                        |              |              3|
| Rajidae            | *Leucoraja fullonica*          | Shagreen ray                     |              |              1|
| Rajidae            | *Neoraja caerulea*             | Blue ray                         |              |              1|
| Rajidae            | *Raja miraletus*               | Brown ray                        |              |              2|
| Rajidae            | *Rajella bathyphila*           | Deep-water ray                   |              |               |
| Rajidae            | *Rajella bigelowi*             | Bigelow's ray                    |              |              1|
| Rajidae            | *Rajella fyllae*               | Round ray                        |              |              2|
| Rajidae            | *Rostroraja alba*              | White skate                      |              |              1|
| Centrophoridae     | *Centrophorus granulosus*      | Gulper shark                     |              |              2|
| Centrophoridae     | *Centrophorus squamosus*       | Leafscale gulper shark           |              |              3|
| Centrophoridae     | *Deania calcea*                | Birdbeak dogfish                 |              |              2|
| Dalatiidae         | *Dalatias licha*               | Kitefin shark                    |              |              4|
| Echinorhinidae     | *Echinorhinus brucus*          | Bramble shark                    |              |               |
| Etmopteridae       | *Centroscyllium fabricii*      | Black dogfish                    |              |              2|
| Etmopteridae       | *Etmopterus princeps*          | Great lanternshark               |              |              1|
| Etmopteridae       | *Etmopterus spinax*            | Velvet belly                     |              |              3|
| Oxynotidae         | *Oxynotus centrina*            | Angular roughshark               |              |               |
| Oxynotidae         | *Oxynotus paradoxus*           | Sailfin roughshark               |              |              1|
| Somniosidae        | *Centroscymnus coelolepis*     | Portuguese dogfish               |              |               |
| Somniosidae        | *Centroscymnus crepidater*     | Longnose velvet dogfish          |              |               |
| Somniosidae        | *Centroscymnus owstonii*       | Roughskin dogfish                |              |               |
| Somniosidae        | *Scymnodon ringens*            | Knifetooth dogfish               |              |              1|
| Somniosidae        | *Somniosus microcephalus*      | Greenland shark                  |              |              4|
| Somniosidae        | *Somniosus rostratus*          | Little sleeper shark             |              |               |
| Somniosidae        | *Zameus squamulosus*           | Velvet dogfish                   |              |               |
| Squatinidae        | *Squatina squatina*            | Angelshark                       |              |              2|
| Torpedinidae       | *Torpedo marmorata*            | Marbled electric ray             |              |              1|
| Torpedinidae       | *Torpedo nobiliana*            | Electric ray                     |              |               |
| Chimaeridae        | *Chimaera monstrosa*           | Rabbit fish                      |              |              4|
| Chimaeridae        | *Chimaera opalescens*          | Opal chimaera                    |              |              4|
| Rhinochimaeridae   | *Harriotta raleighana*         | Pacific longnose chimaera        |              |              3|
| Rhinochimaeridae   | *Rhinochimaera atlantica*      | Straightnose rabbitfish          |              |               |
| Myxinidae          | *Myxine glutinosa*             | Atlantic hagfish                 |              |              2|

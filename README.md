### PSM-Narcolepsy
This repository was designed to assess the relation between Narcolepsy and the Influenza vaccine through the study of protein sequences in two Influenza Virus Vaccine strains and protein sequences in the body.

#### Preparation 
The given data is filtered, taking only the instances of the _HA-_ proteins (19<sup>th<\sup> column of the given data). Only the name of the protein, the sequence, and the initial position are taken

#### Sequence mapping and counting
First, the sequence of the filtered data (thus, only the *HA-* proteins) are mapped for each position between the range set in `options.json` ) "---" ergo, the sequence of Aminoacids (AA) is written for each position between the set position range. Only sequences that overlap at least one position with the defined position range are considered. Furthermore, Posttranslational Modifications (PTM), defined between squared brackets in the sequence (\[\]) have been removed of the sequence (not affecting the indexing of the positions). 

Subsequently, the number of repeated sequences overlapping to some extent with the proposed range of position is counted. Similarly to mapping the sequences, only the AA sequences that overlap at least in one position with the defined range are considered, and PTM was removed. Two outputs are defined: `seqCount` and `seqCountRange`. The propper corresponds to the number of sequences overlapping at least one position with the defined range, considering the positions outside of the range. The latter considers the sequences that overlap at least one position with the range, but only the AAs within the range were considered (thus does not include any AA outside the defined range).

Finally, in order to visually assess the mutations, each sequence is chronologically sorted based on their initial positions, and compared to a reference protein (_refProt_ in `options.json`). The number of repetitions of the sequence for each strain is included at the end of the positions. Mutations have been highlighted

#### Commentary
All the settings of the repository must be set in `options.json`. The storage of the data is of free-will, yet paths must be defined in the `options.json` file.

`mapSequences.py` outputs a _.txt_ file, whereas `mapSequencesHTML.py` outputs an _HTML_ file. `options.json` allows to modify the way the mutations are highlighted under _"MarkdownHigh"_.

### Results 
Results are available in the folder `./Outputs/`, both in _.txt_ format and _.html_. In order to properly visualize the results, please open the proper in a text editor, and the latter in a browser.

         LVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKC

EGRMNYYWTLVEPGDK(ARP:1, PAN:0) // **+15.995**:1

   MNYYWTLVEPGDK(ARP:52, PAN:430) // **+0.984**:36 // **+31.990**:66 // **+15.995**:79 // **+14.016**:10 // **+43.006**:10

   MNYYWTLVEPGDKITFEATGNLVVPR(ARP:28, PAN:48) // **+31.990**:10 // **+15.995**:21 // **+0.984**:5 // **+43.006**:1 // **+28.031**:1

   MNYYWTLVEPGDK<span style="color: red">R</span>(ARP:4, PAN:17) // **+15.995**:8 // **+0.984**:1

   MNYYWT<span style="color: red">Q</span>VEPGDKITFEATGNLVVPR(ARP:2, PAN:1) // **+31.990**:1 // **+0.984**:2 // **+15.995**:2

   MNYYWTLVEPGDK<span style="color: red">R</span>TFEATGNLVVPR(ARP:0, PAN:21) // **+28.031**:21 // **+15.995**:3 // **+0.984**:1

   MNYYWTLVEPGDKITFEATGNLVVPRYAFAMER(ARP:0, PAN:33) // **+15.995**:1 // **+0.984**:1

   MNYYWTLV<span style="color: red">K</span>PGDKITFEATGNLVVPRYAFAMER(ARP:0, PAN:1)

   MNYYWTLV<span style="color: red">K</span>PGDKITFEATGNLVVPR(ARP:0, PAN:2) // **+0.984**:1 // **+79.966**:1

   MNYYWTLV<span style="color: red">K</span>PGDK(ARP:0, PAN:1) // **+31.990**:1

      YWTLVEPGDKITF(ARP:6, PAN:14)

      <span style="color: red">H</span>WTLVEPGDKITF(ARP:0, PAN:1) // **+42.011**:1

      YWT<span style="color: red">Q</span>VEPGDKITF(ARP:0, PAN:1)

      Y<span style="color: red">C</span>TLVEPGDKITF(ARP:0, PAN:2) // **+43.006**:2 // **+71.037**:2

      YWTLVEPGDK<span style="color: red">R</span>TF(ARP:0, PAN:6) // **+28.031**:6

       WTLVEPGDKITF(ARP:8, PAN:16)

       <span style="color: red">R</span>TLVEPGDKITF(ARP:0, PAN:1) // **+43.006**:1

       WT<span style="color: red">P</span>VEPGDKITF(ARP:0, PAN:1) // **+28.031**:1

       W<span style="color: red">A</span>LVEPGDKITF(ARP:0, PAN:1) // **+43.006**:1 // **+42.011**:1

       WTLVEPGDK<span style="color: red">R</span>TF(ARP:0, PAN:6) // **+28.031**:6

       WT<span style="color: red">Q</span>VEPGDKITF(ARP:0, PAN:1)

        TLVEPGDK(ARP:0, PAN:3)

        TLVEPGDKITF(ARP:7, PAN:284) // **+43.006**:1 // **+28.031**:1 // **+14.016**:3

        TLVEPGDK<span style="color: red">R</span>TF(ARP:0, PAN:1) // **+28.031**:1

        TLVEPGDKITFEATGNLVVPRYAF(ARP:3, PAN:6)

        TLVEPGDKITFEATGNLVVPRY(ARP:0, PAN:1)

        TLVEPGDKITF<span style="color: red">Q</span>ATGNLVVPRYAF(ARP:0, PAN:1) // **+0.984**:1

        TLVEPGDK<span style="color: red">V</span>TF(ARP:0, PAN:1) // **+42.011**:1

        TL<span style="color: red">I</span>EPGDKITF(ARP:0, PAN:1)

        TL<span style="color: red">L</span>EPGDKITF(ARP:0, PAN:1)

        TLVEPGDKITFEATG<span style="color: red">Y</span>(ARP:0, PAN:1) // **+43.006**:1

        TLVEPGDKIT<span style="color: red">Y</span>(ARP:0, PAN:4) // **+43.006**:2 // **+42.011**:1 // **+28.031**:1 // **+203.079**:1

        TLVEP<span style="color: red">A</span>DKITF(ARP:0, PAN:1) // **+43.006**:2

        TL<span style="color: red">E</span>EPGDKITF(ARP:0, PAN:3) // **+43.006**:3 // **+28.031**:2 // **+203.079**:1

        TLV<span style="color: red">Q</span>PGDKITF(ARP:0, PAN:1) // **+43.006**:1 // **+79.966**:1

        TL<span style="color: red">A</span>EPGDKITF(ARP:0, PAN:1) // **+43.006**:1 // **+42.011**:1

        TLV<span style="color: red">G</span>PGDKITF(ARP:0, PAN:1) // **+42.011**:1

        TLVEPG<span style="color: red">N</span>KITF(ARP:0, PAN:1) // **+43.006**:1 // **+79.966**:1

        TLVE<span style="color: red">T</span>GDKITF(ARP:0, PAN:1) // **+43.006**:1

        <span style="color: red">I</span>LVEPGDKITF(ARP:0, PAN:1) // **+14.016**:1

        TLVEPG<span style="color: red">A</span>KITF(ARP:0, PAN:1) // **+28.031**:1

        <span style="color: red">K</span>LVEPGDKITF(ARP:0, PAN:1) // **+79.966**:1 // **+43.006**:1

        TLVEPGDKI<span style="color: red">K</span>F(ARP:0, PAN:1) // **+43.006**:1 // **+28.031**:1

        T<span style="color: red">R</span>VEPGDKITF(ARP:0, PAN:1) // **+43.006**:1 // **+28.031**:1

        TLVEPGD<span style="color: red">N</span>ITF(ARP:0, PAN:1) // **+43.006**:1

                ITFEATGNLVVPR(ARP:64, PAN:477) // **+0.984**:21

                ITFEATGNLVVPRYAFAMER(ARP:0, PAN:2)

                <span style="color: red">R</span>TFEATGNLVVPR(ARP:0, PAN:3) // **+43.006**:6

                ITFEATGNLVVPRY<span style="color: red">V</span>FAMER(ARP:0, PAN:3) // **+0.984**:3 // **+43.006**:3

                 TFEATGNLVVPR(ARP:0, PAN:3)

                   EATGNLVVPRY(ARP:4, PAN:16) // **+0.984**:4 // **-18.011**:2

                   EATGNLVVPRYAF(ARP:15, PAN:62) // **+0.984**:14 // **-18.011**:12

                   EATGNLVVPRY<span style="color: red">V</span>F(ARP:5, PAN:1) // **+43.006**:6

                   EATG<span style="color: red">D</span>LVVPRYAF(ARP:0, PAN:3) // **-18.011**:2

                   <span style="color: red">Q</span>ATGNLVVPRYAF(ARP:0, PAN:2) // **+0.984**:2

                   EATG<span style="color: red">K</span>LVVPRY(ARP:0, PAN:1)

                   EAT<span style="color: red">V</span>NLVVPRY(ARP:0, PAN:2) // **+43.006**:2

                   EAT<span style="color: red">E</span>NLVVPRY(ARP:0, PAN:1) // **-18.011**:1 // **+43.006**:1

                   EATGNLVV<span style="color: red">L</span>RY(ARP:0, PAN:1) // **-18.011**:1 // **+0.984**:1

                   EATGNL<span style="color: red">M</span>VPRY(ARP:0, PAN:2)

                   E<span style="color: red">S</span>TGNLVVPRY(ARP:0, PAN:1)

                   EATGNL<span style="color: red">L</span>VPRY(ARP:0, PAN:1) // **+43.006**:2

                   EATGNLVV<span style="color: red">A</span>RY(ARP:0, PAN:1) // **-18.011**:1

                   <span style="color: red">A</span>ATGNLVVPRY(ARP:0, PAN:1) // **+43.006**:1

                   EATGNLVVPR<span style="color: red">F</span>AF(ARP:0, PAN:1)

                   EA<span style="color: red">P</span>GNLVVPRYAF(ARP:0, PAN:1) // **+43.006**:2

                   EA<span style="color: red">A</span>GNLVVPRY(ARP:0, PAN:1) // **+43.006**:1

                   EATGNLVVPRY<span style="color: red">S</span>F(ARP:0, PAN:1) // **-18.011**:1

                   EATGNLV<span style="color: red">E</span>PRYAF(ARP:0, PAN:1)

                   EAT<span style="color: red">R</span>NLVVPRY(ARP:0, PAN:1) // **+43.006**:1 // **+0.984**:1

                   EATGNLVVPR<span style="color: red">D</span>AF(ARP:0, PAN:1) // **-18.011**:1 // **+43.006**:1

                   EATGNLVVP<span style="color: red">S</span>Y(ARP:0, PAN:1)

                        LVVPRYAF(ARP:0, PAN:1)

                             YAFAMER(ARP:52, PAN:113) // **+15.995**:81 // **+31.990**:1

                             YAFAMERNAGSGIIISDTPVHDCNTTCQTPK(ARP:0, PAN:2) // **+71.037**:4 // **+0.984**:1

                             Y<span style="color: red">V</span>FAMER(ARP:0, PAN:9) // **+43.006**:9 // **+15.995**:2

                              AFAMERNAGSGII<span style="color: red">F</span>(ARP:0, PAN:1) // **+43.006**:2

                              AFAMERNAGSG<span style="color: red">F</span>(ARP:0, PAN:1) // **+43.006**:1

                                AMERNAGSGII<span style="color: red">F</span>(ARP:0, PAN:1) // **+43.006**:1 // **+0.984**:1

                                AMERNAGSGIIISDTP<span style="color: red">F</span>(ARP:0, PAN:1)

                                AMERNAGSGIIISDTPVHDCNTTCQTPKGAI<span style="color: red">S</span>TSLPF(ARP:0, PAN:1) // **+71.037**:2 // **+203.079**:2 // **+42.011**:1

                                    NAGSGIIISDTPVHDCNTTCQTPK(ARP:9, PAN:23) // **+71.037**:64 // **+203.079**:17 // **+0.984**:12 // **+28.031**:2

                                           ISDTPVHDCNTTCQTPKGAINTSLPF(ARP:0, PAN:1) // **+71.037**:2 // **+0.984**:1

                                                     TTCQTPKGAINTSLPF(ARP:0, PAN:1) // **+71.037**:1 // **+0.984**:1

                                                        QTPKGAINTSLPF(ARP:0, PAN:1) // **+203.079**:1

                                                            GAINTSLPFQNIHPITIGK(ARP:16, PAN:67) // **+203.079**:75 // **+0.984**:20

                                                            GAINTSLPFQNIH<span style="color: red">Q</span>ITIGK(ARP:13, PAN:4) // **+203.079**:16 // **+0.984**:10

                                                            GAINTSLPFQNIHPITIGKCPK(ARP:1, PAN:1) // **+203.079**:2 // **+71.037**:2 // **+0.984**:1

                                                                     QNIHPITIGKCPKY(ARP:8, PAN:48) // **+71.037**:56 // **-17.027**:28 // **+0.984**:1

                                                                     QNIH<span style="color: red">Q</span>ITIGKCPKY(ARP:5, PAN:0) // **+71.037**:5 // **-17.027**:2

                                                                     QNIHPITIGKC<span style="color: red">A</span>KY(ARP:0, PAN:1) // **+71.037**:1 // **+42.011**:1

                                                                     <span style="color: red">E</span>NIHPITIGKCPKY(ARP:0, PAN:2) // **-18.011**:2 // **+71.037**:2

                                                                     QNI<span style="color: red">D</span>PITIGKCPKY(ARP:0, PAN:1) // **+71.037**:1

                                                                     QNIHPITIGKCP<span style="color: red">T</span>Y(ARP:0, PAN:1) // **-17.027**:1 // **+43.006**:1 // **+71.037**:1

                                                                     QNIHPITI<span style="color: red">R</span>KCPKY(ARP:0, PAN:2) // **+71.037**:2 // **-17.027**:1

                                                                     QNI<span style="color: red">N</span>PITIGKCPKY(ARP:0, PAN:1) // **+43.006**:2 // **+71.037**:1 // **+42.011**:1

                                                                     QNIHPI<span style="color: red">P</span>IGKCPKY(ARP:0, PAN:1) // **+43.006**:2 // **+71.037**:1

                                                                     QNIHPITIGK<span style="color: red">S</span>PKY(ARP:0, PAN:1) // **-17.027**:1 // **+43.006**:2

                                                                     Q<span style="color: red">I</span>IHPITIGKCPKY(ARP:0, PAN:1) // **+79.966**:1 // **+71.037**:1

                                                                     QNIHPITIGKC<span style="color: red">S</span>KY(ARP:0, PAN:1) // **-17.027**:1 // **+43.006**:1 // **+71.037**:1

                                                                     QNIHPI<span style="color: red">A</span>IGKCPKY(ARP:0, PAN:1) // **-17.027**:1 // **+71.037**:1

                                                                         PITIGKCPKY(ARP:0, PAN:1) // **+43.006**:2 // **+71.037**:1
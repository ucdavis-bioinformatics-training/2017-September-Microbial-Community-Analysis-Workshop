Using the Phyloseq package
==========================

installation from bioconductor
------------------------------

``` r
#source("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")
library(phyloseq)
library(ggplot2)
```

Read in the dataset, biom file generated from dbcAmplicons pipeline
-------------------------------------------------------------------

``` r
slashpile_16sV1V3 <- "16sV1V3.biom"
s16sV1V3 = import_biom(BIOMfilename = slashpile_16sV1V3, parseFunction = parse_taxonomy_default)
colnames(tax_table(s16sV1V3)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
rank_names(s16sV1V3)
```

    ## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"

``` r
s16sV1V3
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 950 taxa and 28 samples ]
    ## sample_data() Sample Data:       [ 28 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 950 taxa by 6 taxonomic ranks ]

Filtering
---------

Lets generate a prevelance table for each taxa

``` r
prevelancedf = apply(X = otu_table(s16sV1V3),
                 MARGIN = 1,
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevelancedf = data.frame(Prevalence = prevelancedf,
                      TotalAbundance = taxa_sums(s16sV1V3),
                      tax_table(s16sV1V3))
prevelancedf
```

    ##            Prevalence TotalAbundance     Kingdom
    ## Taxa_00000          8             17  d__Archaea
    ## Taxa_00001          1              1  d__Archaea
    ## Taxa_00002         28          63312 d__Bacteria
    ## Taxa_00003         28           1864 d__Bacteria
    ## Taxa_00004         28           9065 d__Bacteria
    ## Taxa_00005         22            248 d__Bacteria
    ## Taxa_00006          6             45 d__Bacteria
    ## Taxa_00007         19             71 d__Bacteria
    ## Taxa_00008         23            183 d__Bacteria
    ## Taxa_00009         25            352 d__Bacteria
    ## Taxa_00010         28           2791 d__Bacteria
    ## Taxa_00011         25            238 d__Bacteria
    ## Taxa_00012          4              4 d__Bacteria
    ## Taxa_00013          4              4 d__Bacteria
    ## Taxa_00014         25            598 d__Bacteria
    ## Taxa_00015         26            548 d__Bacteria
    ## Taxa_00016          6              7 d__Bacteria
    ## Taxa_00017         27           3029 d__Bacteria
    ## Taxa_00018         24            405 d__Bacteria
    ## Taxa_00019         28            895 d__Bacteria
    ## Taxa_00020         28          26211 d__Bacteria
    ## Taxa_00021         28           8990 d__Bacteria
    ## Taxa_00022         23            541 d__Bacteria
    ## Taxa_00023         22            116 d__Bacteria
    ## Taxa_00024          7              9 d__Bacteria
    ## Taxa_00025         14            151 d__Bacteria
    ## Taxa_00026          9             37 d__Bacteria
    ## Taxa_00027         27          27849 d__Bacteria
    ## Taxa_00028         28           2356 d__Bacteria
    ## Taxa_00029         17             50 d__Bacteria
    ## Taxa_00030         28           3626 d__Bacteria
    ## Taxa_00031         28          14550 d__Bacteria
    ## Taxa_00032         28            347 d__Bacteria
    ## Taxa_00033         27            815 d__Bacteria
    ## Taxa_00034         28           1746 d__Bacteria
    ## Taxa_00035         28            953 d__Bacteria
    ## Taxa_00036         28          11523 d__Bacteria
    ## Taxa_00037          2              2 d__Bacteria
    ## Taxa_00038         27           2110 d__Bacteria
    ## Taxa_00039         28          15373 d__Bacteria
    ## Taxa_00040         28           5550 d__Bacteria
    ## Taxa_00041          2              2 d__Bacteria
    ## Taxa_00042          5             12 d__Bacteria
    ## Taxa_00043          1              1 d__Bacteria
    ## Taxa_00044         27            724 d__Bacteria
    ## Taxa_00045         28           2505 d__Bacteria
    ## Taxa_00046         28           3398 d__Bacteria
    ## Taxa_00047         19             67 d__Bacteria
    ## Taxa_00048          1              1 d__Bacteria
    ## Taxa_00049          3              3 d__Bacteria
    ## Taxa_00050         22            558 d__Bacteria
    ## Taxa_00051         28           2404 d__Bacteria
    ## Taxa_00052         28            789 d__Bacteria
    ## Taxa_00053         28           1740 d__Bacteria
    ## Taxa_00054         27            170 d__Bacteria
    ## Taxa_00055         28           8699 d__Bacteria
    ## Taxa_00056          4              5 d__Bacteria
    ## Taxa_00057          3              4 d__Bacteria
    ## Taxa_00058          1              1 d__Bacteria
    ## Taxa_00059          1              1 d__Bacteria
    ## Taxa_00060         18            120 d__Bacteria
    ## Taxa_00061          1              1 d__Bacteria
    ## Taxa_00062          1              1 d__Bacteria
    ## Taxa_00063          2              2 d__Bacteria
    ## Taxa_00064          8             24 d__Bacteria
    ## Taxa_00065         23             95 d__Bacteria
    ## Taxa_00066         19             63 d__Bacteria
    ## Taxa_00067          5              7 d__Bacteria
    ## Taxa_00068         28            554 d__Bacteria
    ## Taxa_00069          2              3 d__Bacteria
    ## Taxa_00070          2              3 d__Bacteria
    ## Taxa_00071          4              5 d__Bacteria
    ## Taxa_00072          1              1 d__Bacteria
    ## Taxa_00073          7             10 d__Bacteria
    ## Taxa_00074          2              2 d__Bacteria
    ## Taxa_00075          9             37 d__Bacteria
    ## Taxa_00076         14            159 d__Bacteria
    ## Taxa_00077         13             85 d__Bacteria
    ## Taxa_00078         14             40 d__Bacteria
    ## Taxa_00079          1              1 d__Bacteria
    ## Taxa_00080          6            187 d__Bacteria
    ## Taxa_00081          1              1 d__Bacteria
    ## Taxa_00082          6              8 d__Bacteria
    ## Taxa_00083          7             18 d__Bacteria
    ## Taxa_00084          3              6 d__Bacteria
    ## Taxa_00085          1              2 d__Bacteria
    ## Taxa_00086         18            125 d__Bacteria
    ## Taxa_00087         27           1250 d__Bacteria
    ## Taxa_00088          6             21 d__Bacteria
    ## Taxa_00089          2              3 d__Bacteria
    ## Taxa_00090          4              7 d__Bacteria
    ## Taxa_00091          4              6 d__Bacteria
    ## Taxa_00092         11             29 d__Bacteria
    ## Taxa_00093         12            153 d__Bacteria
    ## Taxa_00094          1              1 d__Bacteria
    ## Taxa_00095         11             26 d__Bacteria
    ## Taxa_00096         10             18 d__Bacteria
    ## Taxa_00097          1              1 d__Bacteria
    ## Taxa_00098          9             26 d__Bacteria
    ## Taxa_00099          4             14 d__Bacteria
    ## Taxa_00100          9             29 d__Bacteria
    ## Taxa_00101          2              4 d__Bacteria
    ## Taxa_00102          2              3 d__Bacteria
    ## Taxa_00103          2              2 d__Bacteria
    ## Taxa_00104         14            101 d__Bacteria
    ## Taxa_00105         15            116 d__Bacteria
    ## Taxa_00106          1              2 d__Bacteria
    ## Taxa_00107          1             10 d__Bacteria
    ## Taxa_00108          7            173 d__Bacteria
    ## Taxa_00109          5             20 d__Bacteria
    ## Taxa_00110         10            141 d__Bacteria
    ## Taxa_00111          1              1 d__Bacteria
    ## Taxa_00112         11             38 d__Bacteria
    ## Taxa_00113          5              6 d__Bacteria
    ## Taxa_00114          4              7 d__Bacteria
    ## Taxa_00115          1              3 d__Bacteria
    ## Taxa_00116         21            809 d__Bacteria
    ## Taxa_00117          2              4 d__Bacteria
    ## Taxa_00118         25           4343 d__Bacteria
    ## Taxa_00119          3              3 d__Bacteria
    ## Taxa_00120          6             17 d__Bacteria
    ## Taxa_00121          1              1 d__Bacteria
    ## Taxa_00122         28           1206 d__Bacteria
    ## Taxa_00123         17             97 d__Bacteria
    ## Taxa_00124          3              3 d__Bacteria
    ## Taxa_00125          4              5 d__Bacteria
    ## Taxa_00126          1              1 d__Bacteria
    ## Taxa_00127          1              1 d__Bacteria
    ## Taxa_00128          1              1 d__Bacteria
    ## Taxa_00129         20             75 d__Bacteria
    ## Taxa_00130          4              4 d__Bacteria
    ## Taxa_00131          9             16 d__Bacteria
    ## Taxa_00132          9             21 d__Bacteria
    ## Taxa_00133          4              5 d__Bacteria
    ## Taxa_00134          5             10 d__Bacteria
    ## Taxa_00135         11             20 d__Bacteria
    ## Taxa_00136         13             23 d__Bacteria
    ## Taxa_00137          4              7 d__Bacteria
    ## Taxa_00138         13            339 d__Bacteria
    ## Taxa_00139          5              8 d__Bacteria
    ## Taxa_00140          4              5 d__Bacteria
    ## Taxa_00141         12             28 d__Bacteria
    ## Taxa_00142          1              1 d__Bacteria
    ## Taxa_00143          3              9 d__Bacteria
    ## Taxa_00144         16             86 d__Bacteria
    ## Taxa_00145          7             23 d__Bacteria
    ## Taxa_00146         11             29 d__Bacteria
    ## Taxa_00147          9             14 d__Bacteria
    ## Taxa_00148         28           3493 d__Bacteria
    ## Taxa_00149         21            175 d__Bacteria
    ## Taxa_00150         22            253 d__Bacteria
    ## Taxa_00151         10             14 d__Bacteria
    ## Taxa_00152         22            319 d__Bacteria
    ## Taxa_00153          6            103 d__Bacteria
    ## Taxa_00154         23            522 d__Bacteria
    ## Taxa_00155          5             10 d__Bacteria
    ## Taxa_00156          1              1 d__Bacteria
    ## Taxa_00157         11             35 d__Bacteria
    ## Taxa_00158         10             40 d__Bacteria
    ## Taxa_00159         15            104 d__Bacteria
    ## Taxa_00160         17            115 d__Bacteria
    ## Taxa_00161         11             17 d__Bacteria
    ## Taxa_00162         11            132 d__Bacteria
    ## Taxa_00163          3              8 d__Bacteria
    ## Taxa_00164         12             17 d__Bacteria
    ## Taxa_00165          1              1 d__Bacteria
    ## Taxa_00166          1              1 d__Bacteria
    ## Taxa_00167         15             37 d__Bacteria
    ## Taxa_00168          1              1 d__Bacteria
    ## Taxa_00169          6             14 d__Bacteria
    ## Taxa_00170          7             22 d__Bacteria
    ## Taxa_00171          2              2 d__Bacteria
    ## Taxa_00172          4              7 d__Bacteria
    ## Taxa_00173         24            302 d__Bacteria
    ## Taxa_00174          1              1 d__Bacteria
    ## Taxa_00175          2              2 d__Bacteria
    ## Taxa_00176          3              3 d__Bacteria
    ## Taxa_00177          1              1 d__Bacteria
    ## Taxa_00178          1              1 d__Bacteria
    ## Taxa_00179         25            349 d__Bacteria
    ## Taxa_00180          3              8 d__Bacteria
    ## Taxa_00181          2              3 d__Bacteria
    ## Taxa_00182         12             23 d__Bacteria
    ## Taxa_00183          8             12 d__Bacteria
    ## Taxa_00184          1              1 d__Bacteria
    ## Taxa_00185          2              2 d__Bacteria
    ## Taxa_00186         21            225 d__Bacteria
    ## Taxa_00187          3              4 d__Bacteria
    ## Taxa_00188          1              1 d__Bacteria
    ## Taxa_00189         11             23 d__Bacteria
    ## Taxa_00190          8             11 d__Bacteria
    ## Taxa_00191          1              1 d__Bacteria
    ## Taxa_00192          2              2 d__Bacteria
    ## Taxa_00193         24            118 d__Bacteria
    ## Taxa_00194         26            459 d__Bacteria
    ## Taxa_00195         25            303 d__Bacteria
    ## Taxa_00196         21            195 d__Bacteria
    ## Taxa_00197         27            905 d__Bacteria
    ## Taxa_00198         18            102 d__Bacteria
    ## Taxa_00199          2              2 d__Bacteria
    ## Taxa_00200          3              3 d__Bacteria
    ## Taxa_00201          1              3 d__Bacteria
    ## Taxa_00202          7             13 d__Bacteria
    ## Taxa_00203         12             39 d__Bacteria
    ## Taxa_00204          2              4 d__Bacteria
    ## Taxa_00205          4              4 d__Bacteria
    ## Taxa_00206         20             89 d__Bacteria
    ## Taxa_00207         17            110 d__Bacteria
    ## Taxa_00208          7              7 d__Bacteria
    ## Taxa_00209         20            181 d__Bacteria
    ## Taxa_00210          4              4 d__Bacteria
    ## Taxa_00211          1              3 d__Bacteria
    ## Taxa_00212         28           7250 d__Bacteria
    ## Taxa_00213         28            809 d__Bacteria
    ## Taxa_00214         28            695 d__Bacteria
    ## Taxa_00215         23            105 d__Bacteria
    ## Taxa_00216         28           1837 d__Bacteria
    ## Taxa_00217         28            280 d__Bacteria
    ## Taxa_00218          2              2 d__Bacteria
    ## Taxa_00219          1              1 d__Bacteria
    ## Taxa_00220          1              1 d__Bacteria
    ## Taxa_00221         25             86 d__Bacteria
    ## Taxa_00222          5              8 d__Bacteria
    ## Taxa_00223         26            871 d__Bacteria
    ## Taxa_00224         28            757 d__Bacteria
    ## Taxa_00225         28            386 d__Bacteria
    ## Taxa_00226         28           1839 d__Bacteria
    ## Taxa_00227         27            416 d__Bacteria
    ## Taxa_00228         28           4260 d__Bacteria
    ## Taxa_00229         21            104 d__Bacteria
    ## Taxa_00230          6              7 d__Bacteria
    ## Taxa_00231          3              5 d__Bacteria
    ## Taxa_00232          1              2 d__Bacteria
    ## Taxa_00233          3              3 d__Bacteria
    ## Taxa_00234         27           1456 d__Bacteria
    ## Taxa_00235         28           1796 d__Bacteria
    ## Taxa_00236         10             26 d__Bacteria
    ## Taxa_00237          1              2 d__Bacteria
    ## Taxa_00238          1              1 d__Bacteria
    ## Taxa_00239         19             63 d__Bacteria
    ## Taxa_00240          5           1615 d__Bacteria
    ## Taxa_00241          3              3 d__Bacteria
    ## Taxa_00242          2            105 d__Bacteria
    ## Taxa_00243          1              5 d__Bacteria
    ## Taxa_00244          1              2 d__Bacteria
    ## Taxa_00245          5             12 d__Bacteria
    ## Taxa_00246         20             68 d__Bacteria
    ## Taxa_00247          1              1 d__Bacteria
    ## Taxa_00248          4              6 d__Bacteria
    ## Taxa_00249         23             91 d__Bacteria
    ## Taxa_00250         17             45 d__Bacteria
    ## Taxa_00251          6              8 d__Bacteria
    ## Taxa_00252          1              1 d__Bacteria
    ## Taxa_00253          1              1 d__Bacteria
    ## Taxa_00254          1              1 d__Bacteria
    ## Taxa_00255          1              1 d__Bacteria
    ## Taxa_00256         28           3271 d__Bacteria
    ## Taxa_00257         28            664 d__Bacteria
    ## Taxa_00258         15             23 d__Bacteria
    ## Taxa_00259          3              4 d__Bacteria
    ## Taxa_00260          6             21 d__Bacteria
    ## Taxa_00261          1              1 d__Bacteria
    ## Taxa_00262          1              1 d__Bacteria
    ## Taxa_00263         27            981 d__Bacteria
    ## Taxa_00264          1              1 d__Bacteria
    ## Taxa_00265         17            126 d__Bacteria
    ## Taxa_00266          1              1 d__Bacteria
    ## Taxa_00267          3              5 d__Bacteria
    ## Taxa_00268          1              1 d__Bacteria
    ## Taxa_00269          1              1 d__Bacteria
    ## Taxa_00270         28            285 d__Bacteria
    ## Taxa_00271         28           6835 d__Bacteria
    ## Taxa_00272          2              2 d__Bacteria
    ## Taxa_00273          9             22 d__Bacteria
    ## Taxa_00274         22            236 d__Bacteria
    ## Taxa_00275         28           1017 d__Bacteria
    ## Taxa_00276          1              1 d__Bacteria
    ## Taxa_00277         28           2225 d__Bacteria
    ## Taxa_00278         24            196 d__Bacteria
    ## Taxa_00279          9            103 d__Bacteria
    ## Taxa_00280         25            273 d__Bacteria
    ## Taxa_00281         28           2287 d__Bacteria
    ## Taxa_00282         22             74 d__Bacteria
    ## Taxa_00283          5              5 d__Bacteria
    ## Taxa_00284          1              1 d__Bacteria
    ## Taxa_00285          4             18 d__Bacteria
    ## Taxa_00286         23            165 d__Bacteria
    ## Taxa_00287         23             93 d__Bacteria
    ## Taxa_00288         17             63 d__Bacteria
    ## Taxa_00289         24            459 d__Bacteria
    ## Taxa_00290         25            465 d__Bacteria
    ## Taxa_00291          2              2 d__Bacteria
    ## Taxa_00292         28           1359 d__Bacteria
    ## Taxa_00293         16             52 d__Bacteria
    ## Taxa_00294         10             20 d__Bacteria
    ## Taxa_00295          2              2 d__Bacteria
    ## Taxa_00296          3             14 d__Bacteria
    ## Taxa_00297         24            508 d__Bacteria
    ## Taxa_00298         12             78 d__Bacteria
    ## Taxa_00299         28           3479 d__Bacteria
    ## Taxa_00300         17             35 d__Bacteria
    ## Taxa_00301         12            596 d__Bacteria
    ## Taxa_00302         13             26 d__Bacteria
    ## Taxa_00303          9             20 d__Bacteria
    ## Taxa_00304         26             95 d__Bacteria
    ## Taxa_00305         28           6988 d__Bacteria
    ## Taxa_00306         27           1083 d__Bacteria
    ## Taxa_00307          1              1 d__Bacteria
    ## Taxa_00308         28           3335 d__Bacteria
    ## Taxa_00309          5              5 d__Bacteria
    ## Taxa_00310          9             18 d__Bacteria
    ## Taxa_00311         11             27 d__Bacteria
    ## Taxa_00312         13             20 d__Bacteria
    ## Taxa_00313          6             11 d__Bacteria
    ## Taxa_00314         27           1853 d__Bacteria
    ## Taxa_00315         26            719 d__Bacteria
    ## Taxa_00316         13             36 d__Bacteria
    ## Taxa_00317          5             10 d__Bacteria
    ## Taxa_00318          1              1 d__Bacteria
    ## Taxa_00319         15             57 d__Bacteria
    ## Taxa_00320          5              8 d__Bacteria
    ## Taxa_00321          6              6 d__Bacteria
    ## Taxa_00322         10             20 d__Bacteria
    ## Taxa_00323         12             33 d__Bacteria
    ## Taxa_00324         16             65 d__Bacteria
    ## Taxa_00325         24            444 d__Bacteria
    ## Taxa_00326          3              4 d__Bacteria
    ## Taxa_00327          7             26 d__Bacteria
    ## Taxa_00328          5              7 d__Bacteria
    ## Taxa_00329          2              6 d__Bacteria
    ## Taxa_00330         21             67 d__Bacteria
    ## Taxa_00331         21             74 d__Bacteria
    ## Taxa_00332          5              7 d__Bacteria
    ## Taxa_00333         25            458 d__Bacteria
    ## Taxa_00334         26           1139 d__Bacteria
    ## Taxa_00335         25           1655 d__Bacteria
    ## Taxa_00336          1              1 d__Bacteria
    ## Taxa_00337         21             73 d__Bacteria
    ## Taxa_00338         18             42 d__Bacteria
    ## Taxa_00339         22             69 d__Bacteria
    ## Taxa_00340          1              1 d__Bacteria
    ## Taxa_00341          1              1 d__Bacteria
    ## Taxa_00342         21             75 d__Bacteria
    ## Taxa_00343          6             32 d__Bacteria
    ## Taxa_00344          4              4 d__Bacteria
    ## Taxa_00345          1              1 d__Bacteria
    ## Taxa_00346          8             16 d__Bacteria
    ## Taxa_00347         26            481 d__Bacteria
    ## Taxa_00348          2              3 d__Bacteria
    ## Taxa_00349          3              4 d__Bacteria
    ## Taxa_00350         15             32 d__Bacteria
    ## Taxa_00351          7             12 d__Bacteria
    ## Taxa_00352          6              6 d__Bacteria
    ## Taxa_00353         28           4442 d__Bacteria
    ## Taxa_00354         27            505 d__Bacteria
    ## Taxa_00355         28           6997 d__Bacteria
    ## Taxa_00356         25            764 d__Bacteria
    ## Taxa_00357         28          17962 d__Bacteria
    ## Taxa_00358         13            892 d__Bacteria
    ## Taxa_00359          3              7 d__Bacteria
    ## Taxa_00360         23            492 d__Bacteria
    ## Taxa_00361         26            664 d__Bacteria
    ## Taxa_00362          1              1 d__Bacteria
    ## Taxa_00363         27           6980 d__Bacteria
    ## Taxa_00364         13            114 d__Bacteria
    ## Taxa_00365         10            241 d__Bacteria
    ## Taxa_00366          3              4 d__Bacteria
    ## Taxa_00367          3              5 d__Bacteria
    ## Taxa_00368         26           1213 d__Bacteria
    ## Taxa_00369          4              4 d__Bacteria
    ## Taxa_00370          1              1 d__Bacteria
    ## Taxa_00371         13             45 d__Bacteria
    ## Taxa_00372          3              3 d__Bacteria
    ## Taxa_00373          1              1 d__Bacteria
    ## Taxa_00374          1              1 d__Bacteria
    ## Taxa_00375          4            161 d__Bacteria
    ## Taxa_00376          9             34 d__Bacteria
    ## Taxa_00377         18            729 d__Bacteria
    ## Taxa_00378         11             35 d__Bacteria
    ## Taxa_00379          1              1 d__Bacteria
    ## Taxa_00380          3              4 d__Bacteria
    ## Taxa_00381          3              3 d__Bacteria
    ## Taxa_00382          1              1 d__Bacteria
    ## Taxa_00383          5             11 d__Bacteria
    ## Taxa_00384         28           1683 d__Bacteria
    ## Taxa_00385          6             14 d__Bacteria
    ## Taxa_00386         24            289 d__Bacteria
    ## Taxa_00387         28           1490 d__Bacteria
    ## Taxa_00388         28           1839 d__Bacteria
    ## Taxa_00389          3              3 d__Bacteria
    ## Taxa_00390          1              1 d__Bacteria
    ## Taxa_00391          1              2 d__Bacteria
    ## Taxa_00392          3             70 d__Bacteria
    ## Taxa_00393          5             10 d__Bacteria
    ## Taxa_00394          7             30 d__Bacteria
    ## Taxa_00395         17           1922 d__Bacteria
    ## Taxa_00396          8            309 d__Bacteria
    ## Taxa_00397          2             49 d__Bacteria
    ## Taxa_00398         16            329 d__Bacteria
    ## Taxa_00399          3             41 d__Bacteria
    ## Taxa_00400          9            596 d__Bacteria
    ## Taxa_00401          9             54 d__Bacteria
    ## Taxa_00402         18           1989 d__Bacteria
    ## Taxa_00403         15            920 d__Bacteria
    ## Taxa_00404          1             16 d__Bacteria
    ## Taxa_00405          5            344 d__Bacteria
    ## Taxa_00406         18           1486 d__Bacteria
    ## Taxa_00407          3              3 d__Bacteria
    ## Taxa_00408          1              2 d__Bacteria
    ## Taxa_00409         19            303 d__Bacteria
    ## Taxa_00410          2             12 d__Bacteria
    ## Taxa_00411          1              1 d__Bacteria
    ## Taxa_00412          8             29 d__Bacteria
    ## Taxa_00413          5              6 d__Bacteria
    ## Taxa_00414         15             31 d__Bacteria
    ## Taxa_00415          8            236 d__Bacteria
    ## Taxa_00416          4              8 d__Bacteria
    ## Taxa_00417          9            747 d__Bacteria
    ## Taxa_00418         12             20 d__Bacteria
    ## Taxa_00419          1              3 d__Bacteria
    ## Taxa_00420          1              1 d__Bacteria
    ## Taxa_00421         12             18 d__Bacteria
    ## Taxa_00422          3              4 d__Bacteria
    ## Taxa_00423          5              6 d__Bacteria
    ## Taxa_00424         27            201 d__Bacteria
    ## Taxa_00425         24            255 d__Bacteria
    ## Taxa_00426          1              2 d__Bacteria
    ## Taxa_00427         16            126 d__Bacteria
    ## Taxa_00428          5              7 d__Bacteria
    ## Taxa_00429         17            167 d__Bacteria
    ## Taxa_00430          7             37 d__Bacteria
    ## Taxa_00431          1              1 d__Bacteria
    ## Taxa_00432          3              9 d__Bacteria
    ## Taxa_00433         12             37 d__Bacteria
    ## Taxa_00434          1              1 d__Bacteria
    ## Taxa_00435          1              1 d__Bacteria
    ## Taxa_00436          1             39 d__Bacteria
    ## Taxa_00437          2              3 d__Bacteria
    ## Taxa_00438          3              3 d__Bacteria
    ## Taxa_00439         24            521 d__Bacteria
    ## Taxa_00440          2              2 d__Bacteria
    ## Taxa_00441          1              3 d__Bacteria
    ## Taxa_00442          4             10 d__Bacteria
    ## Taxa_00443          8             31 d__Bacteria
    ## Taxa_00444          4             10 d__Bacteria
    ## Taxa_00445          1              1 d__Bacteria
    ## Taxa_00446          3              4 d__Bacteria
    ## Taxa_00447          2              2 d__Bacteria
    ## Taxa_00448         14             34 d__Bacteria
    ## Taxa_00449          9             52 d__Bacteria
    ## Taxa_00450          2              2 d__Bacteria
    ## Taxa_00451          3              4 d__Bacteria
    ## Taxa_00452          2              3 d__Bacteria
    ## Taxa_00453          4              5 d__Bacteria
    ## Taxa_00454          5             13 d__Bacteria
    ## Taxa_00455          1              1 d__Bacteria
    ## Taxa_00456          1              2 d__Bacteria
    ## Taxa_00457          6             10 d__Bacteria
    ## Taxa_00458          3              5 d__Bacteria
    ## Taxa_00459          4             92 d__Bacteria
    ## Taxa_00460          7             11 d__Bacteria
    ## Taxa_00461          5              5 d__Bacteria
    ## Taxa_00462          7              8 d__Bacteria
    ## Taxa_00463          2              3 d__Bacteria
    ## Taxa_00464          8             22 d__Bacteria
    ## Taxa_00465          8             53 d__Bacteria
    ## Taxa_00466         16            137 d__Bacteria
    ## Taxa_00467         28          22581 d__Bacteria
    ## Taxa_00468          9             16 d__Bacteria
    ## Taxa_00469         25            398 d__Bacteria
    ## Taxa_00470          4             15 d__Bacteria
    ## Taxa_00471          6             12 d__Bacteria
    ## Taxa_00472         21            192 d__Bacteria
    ## Taxa_00473         22           1387 d__Bacteria
    ## Taxa_00474          3              3 d__Bacteria
    ## Taxa_00475         28           5848 d__Bacteria
    ## Taxa_00476         28            176 d__Bacteria
    ## Taxa_00477         20             34 d__Bacteria
    ## Taxa_00478         27            482 d__Bacteria
    ## Taxa_00479         26            272 d__Bacteria
    ## Taxa_00480         24            204 d__Bacteria
    ## Taxa_00481         13             24 d__Bacteria
    ## Taxa_00482          7              8 d__Bacteria
    ## Taxa_00483         28           9051 d__Bacteria
    ## Taxa_00484         28           2119 d__Bacteria
    ## Taxa_00485         14             29 d__Bacteria
    ## Taxa_00486         28           1752 d__Bacteria
    ## Taxa_00487         28           1385 d__Bacteria
    ## Taxa_00488         12             27 d__Bacteria
    ## Taxa_00489         27           1337 d__Bacteria
    ## Taxa_00490         26            106 d__Bacteria
    ## Taxa_00491          4              9 d__Bacteria
    ## Taxa_00492          5              6 d__Bacteria
    ## Taxa_00493         28           1813 d__Bacteria
    ## Taxa_00494         28           3436 d__Bacteria
    ## Taxa_00495         28           2201 d__Bacteria
    ## Taxa_00496         13             48 d__Bacteria
    ## Taxa_00497         28            775 d__Bacteria
    ## Taxa_00498         28          16013 d__Bacteria
    ## Taxa_00499         28          15381 d__Bacteria
    ## Taxa_00500         16             29 d__Bacteria
    ## Taxa_00501          8             11 d__Bacteria
    ## Taxa_00502          3              5 d__Bacteria
    ## Taxa_00503         28           3797 d__Bacteria
    ## Taxa_00504         11             29 d__Bacteria
    ## Taxa_00505         22            171 d__Bacteria
    ## Taxa_00506         16             58 d__Bacteria
    ## Taxa_00507         22           1071 d__Bacteria
    ## Taxa_00508         28           1106 d__Bacteria
    ## Taxa_00509         28           1574 d__Bacteria
    ## Taxa_00510          2              3 d__Bacteria
    ## Taxa_00511          1              1 d__Bacteria
    ## Taxa_00512          1              1 d__Bacteria
    ## Taxa_00513          2              2 d__Bacteria
    ## Taxa_00514         27            304 d__Bacteria
    ## Taxa_00515         28          34954 d__Bacteria
    ## Taxa_00516          2              2 d__Bacteria
    ## Taxa_00517         28           2535 d__Bacteria
    ## Taxa_00518         20            457 d__Bacteria
    ## Taxa_00519          6             14 d__Bacteria
    ## Taxa_00520         16             57 d__Bacteria
    ## Taxa_00521         22           1555 d__Bacteria
    ## Taxa_00522         23             83 d__Bacteria
    ## Taxa_00523         26            395 d__Bacteria
    ## Taxa_00524         12             47 d__Bacteria
    ## Taxa_00525         24            187 d__Bacteria
    ## Taxa_00526          1              1 d__Bacteria
    ## Taxa_00527         28           9729 d__Bacteria
    ## Taxa_00528         28          26786 d__Bacteria
    ## Taxa_00529          6             50 d__Bacteria
    ## Taxa_00530         28          21120 d__Bacteria
    ## Taxa_00531         28           1764 d__Bacteria
    ## Taxa_00532         28            429 d__Bacteria
    ## Taxa_00533         10             27 d__Bacteria
    ## Taxa_00534         20            213 d__Bacteria
    ## Taxa_00535          1              1 d__Bacteria
    ## Taxa_00536         28           1259 d__Bacteria
    ## Taxa_00537          1              3 d__Bacteria
    ## Taxa_00538          1              4 d__Bacteria
    ## Taxa_00539          1              1 d__Bacteria
    ## Taxa_00540          1              5 d__Bacteria
    ## Taxa_00541         28           1578 d__Bacteria
    ## Taxa_00542          1              1 d__Bacteria
    ## Taxa_00543          3              3 d__Bacteria
    ## Taxa_00544         25            230 d__Bacteria
    ## Taxa_00545          6             10 d__Bacteria
    ## Taxa_00546         25            808 d__Bacteria
    ## Taxa_00547         19             58 d__Bacteria
    ## Taxa_00548         27           1409 d__Bacteria
    ## Taxa_00549         28            772 d__Bacteria
    ## Taxa_00550         17            111 d__Bacteria
    ## Taxa_00551         20            575 d__Bacteria
    ## Taxa_00552         28           3304 d__Bacteria
    ## Taxa_00553         10             16 d__Bacteria
    ## Taxa_00554         17             72 d__Bacteria
    ## Taxa_00555          9             11 d__Bacteria
    ## Taxa_00556         27            405 d__Bacteria
    ## Taxa_00557         12             14 d__Bacteria
    ## Taxa_00558         23            145 d__Bacteria
    ## Taxa_00559          1              1 d__Bacteria
    ## Taxa_00560         23             92 d__Bacteria
    ## Taxa_00561          3              3 d__Bacteria
    ## Taxa_00562         27            205 d__Bacteria
    ## Taxa_00563          2              3 d__Bacteria
    ## Taxa_00564         25            381 d__Bacteria
    ## Taxa_00565          2             38 d__Bacteria
    ## Taxa_00566          1             11 d__Bacteria
    ## Taxa_00567          8              8 d__Bacteria
    ## Taxa_00568          8             85 d__Bacteria
    ## Taxa_00569          3              4 d__Bacteria
    ## Taxa_00570         28           1461 d__Bacteria
    ## Taxa_00571          5             15 d__Bacteria
    ## Taxa_00572         19             78 d__Bacteria
    ## Taxa_00573         24            170 d__Bacteria
    ## Taxa_00574         21            110 d__Bacteria
    ## Taxa_00575          9             43 d__Bacteria
    ## Taxa_00576          5             58 d__Bacteria
    ## Taxa_00577         21            209 d__Bacteria
    ## Taxa_00578         26            451 d__Bacteria
    ## Taxa_00579          6             26 d__Bacteria
    ## Taxa_00580         28           1925 d__Bacteria
    ## Taxa_00581         28           2761 d__Bacteria
    ## Taxa_00582          3              3 d__Bacteria
    ## Taxa_00583         20             81 d__Bacteria
    ## Taxa_00584         28           2147 d__Bacteria
    ## Taxa_00585         13            161 d__Bacteria
    ## Taxa_00586         25            299 d__Bacteria
    ## Taxa_00587         13             41 d__Bacteria
    ## Taxa_00588          4              5 d__Bacteria
    ## Taxa_00589          1              1 d__Bacteria
    ## Taxa_00590         23            768 d__Bacteria
    ## Taxa_00591         17            154 d__Bacteria
    ## Taxa_00592         19             85 d__Bacteria
    ## Taxa_00593          1              1 d__Bacteria
    ## Taxa_00594          1              1 d__Bacteria
    ## Taxa_00595         24            720 d__Bacteria
    ## Taxa_00596          1              1 d__Bacteria
    ## Taxa_00597         23            193 d__Bacteria
    ## Taxa_00598         27            663 d__Bacteria
    ## Taxa_00599         22            189 d__Bacteria
    ## Taxa_00600         13             48 d__Bacteria
    ## Taxa_00601          1              1 d__Bacteria
    ## Taxa_00602          2              2 d__Bacteria
    ## Taxa_00603         12             22 d__Bacteria
    ## Taxa_00604          3              8 d__Bacteria
    ## Taxa_00605          4              4 d__Bacteria
    ## Taxa_00606          6              9 d__Bacteria
    ## Taxa_00607         28           3934 d__Bacteria
    ## Taxa_00608         28           6548 d__Bacteria
    ## Taxa_00609          1              1 d__Bacteria
    ## Taxa_00610          5              7 d__Bacteria
    ## Taxa_00611         27            282 d__Bacteria
    ## Taxa_00612          8             31 d__Bacteria
    ## Taxa_00613          2              2 d__Bacteria
    ## Taxa_00614          1              1 d__Bacteria
    ## Taxa_00615          1              2 d__Bacteria
    ## Taxa_00616          2              2 d__Bacteria
    ## Taxa_00617         18             33 d__Bacteria
    ## Taxa_00618          2              2 d__Bacteria
    ## Taxa_00619         11             22 d__Bacteria
    ## Taxa_00620          4             13 d__Bacteria
    ## Taxa_00621         11             26 d__Bacteria
    ## Taxa_00622          2              2 d__Bacteria
    ## Taxa_00623          1              1 d__Bacteria
    ## Taxa_00624          3              3 d__Bacteria
    ## Taxa_00625         28           1236 d__Bacteria
    ## Taxa_00626          6             13 d__Bacteria
    ## Taxa_00627          6             17 d__Bacteria
    ## Taxa_00628         12             91 d__Bacteria
    ## Taxa_00629         22            100 d__Bacteria
    ## Taxa_00630          1              1 d__Bacteria
    ## Taxa_00631         13             33 d__Bacteria
    ## Taxa_00632          2              2 d__Bacteria
    ## Taxa_00633         28           5665 d__Bacteria
    ## Taxa_00634         28           3839 d__Bacteria
    ## Taxa_00635          5             19 d__Bacteria
    ## Taxa_00636          9             13 d__Bacteria
    ## Taxa_00637          1              1 d__Bacteria
    ## Taxa_00638         28           1130 d__Bacteria
    ## Taxa_00639          7             10 d__Bacteria
    ## Taxa_00640          1              1 d__Bacteria
    ## Taxa_00641          4              5 d__Bacteria
    ## Taxa_00642          1              1 d__Bacteria
    ## Taxa_00643          1              1 d__Bacteria
    ## Taxa_00644          2              2 d__Bacteria
    ## Taxa_00645         27            401 d__Bacteria
    ## Taxa_00646          6              9 d__Bacteria
    ## Taxa_00647          2              2 d__Bacteria
    ## Taxa_00648         26           1502 d__Bacteria
    ## Taxa_00649          2              6 d__Bacteria
    ## Taxa_00650         25            274 d__Bacteria
    ## Taxa_00651         13             32 d__Bacteria
    ## Taxa_00652          1              1 d__Bacteria
    ## Taxa_00653         14             24 d__Bacteria
    ## Taxa_00654          4             10 d__Bacteria
    ## Taxa_00655         11             14 d__Bacteria
    ## Taxa_00656         19             54 d__Bacteria
    ## Taxa_00657          6              8 d__Bacteria
    ## Taxa_00658          7             17 d__Bacteria
    ## Taxa_00659          2              2 d__Bacteria
    ## Taxa_00660         19            113 d__Bacteria
    ## Taxa_00661          3              3 d__Bacteria
    ## Taxa_00662          3             16 d__Bacteria
    ## Taxa_00663          4              4 d__Bacteria
    ## Taxa_00664          1              1 d__Bacteria
    ## Taxa_00665         26           3456 d__Bacteria
    ## Taxa_00666         10            125 d__Bacteria
    ## Taxa_00667         22            138 d__Bacteria
    ## Taxa_00668         21            252 d__Bacteria
    ## Taxa_00669          3              6 d__Bacteria
    ## Taxa_00670         11             46 d__Bacteria
    ## Taxa_00671         16             70 d__Bacteria
    ## Taxa_00672         16             86 d__Bacteria
    ## Taxa_00673         23            320 d__Bacteria
    ## Taxa_00674          8             16 d__Bacteria
    ## Taxa_00675          5              7 d__Bacteria
    ## Taxa_00676          8            244 d__Bacteria
    ## Taxa_00677          2              2 d__Bacteria
    ## Taxa_00678         26           1356 d__Bacteria
    ## Taxa_00679         11             68 d__Bacteria
    ## Taxa_00680          1              2 d__Bacteria
    ## Taxa_00681         18            464 d__Bacteria
    ## Taxa_00682         28          20808 d__Bacteria
    ## Taxa_00683         28           8011 d__Bacteria
    ## Taxa_00684         11             78 d__Bacteria
    ## Taxa_00685          1              1 d__Bacteria
    ## Taxa_00686          1              4 d__Bacteria
    ## Taxa_00687          2             18 d__Bacteria
    ## Taxa_00688          5              9 d__Bacteria
    ## Taxa_00689          1             42 d__Bacteria
    ## Taxa_00690          1              1 d__Bacteria
    ## Taxa_00691          1              3 d__Bacteria
    ## Taxa_00692          4              5 d__Bacteria
    ## Taxa_00693          1             36 d__Bacteria
    ## Taxa_00694         19            186 d__Bacteria
    ## Taxa_00695         24           2105 d__Bacteria
    ## Taxa_00696          6             10 d__Bacteria
    ## Taxa_00697         17            117 d__Bacteria
    ## Taxa_00698          2              2 d__Bacteria
    ## Taxa_00699         10             43 d__Bacteria
    ## Taxa_00700         28           1540 d__Bacteria
    ## Taxa_00701          4              4 d__Bacteria
    ## Taxa_00702          6              8 d__Bacteria
    ## Taxa_00703         25            227 d__Bacteria
    ## Taxa_00704          5              8 d__Bacteria
    ## Taxa_00705         27           1420 d__Bacteria
    ## Taxa_00706         11             12 d__Bacteria
    ## Taxa_00707         28            180 d__Bacteria
    ## Taxa_00708         10             11 d__Bacteria
    ## Taxa_00709          2              4 d__Bacteria
    ## Taxa_00710         27            276 d__Bacteria
    ## Taxa_00711          4              5 d__Bacteria
    ## Taxa_00712         28            671 d__Bacteria
    ## Taxa_00713          5              5 d__Bacteria
    ## Taxa_00714          1              1 d__Bacteria
    ## Taxa_00715         19            349 d__Bacteria
    ## Taxa_00716          1              1 d__Bacteria
    ## Taxa_00717         17             98 d__Bacteria
    ## Taxa_00718          2              2 d__Bacteria
    ## Taxa_00719          1              4 d__Bacteria
    ## Taxa_00720          5              6 d__Bacteria
    ## Taxa_00721          2              2 d__Bacteria
    ## Taxa_00722          1              1 d__Bacteria
    ## Taxa_00723         19             52 d__Bacteria
    ## Taxa_00724          4             11 d__Bacteria
    ## Taxa_00725         17            267 d__Bacteria
    ## Taxa_00726          2              2 d__Bacteria
    ## Taxa_00727          2              5 d__Bacteria
    ## Taxa_00728         26            770 d__Bacteria
    ## Taxa_00729         21            164 d__Bacteria
    ## Taxa_00730         10             19 d__Bacteria
    ## Taxa_00731          1              1 d__Bacteria
    ## Taxa_00732         14            277 d__Bacteria
    ## Taxa_00733          4              4 d__Bacteria
    ## Taxa_00734         15             24 d__Bacteria
    ## Taxa_00735         28           4408 d__Bacteria
    ## Taxa_00736         21            219 d__Bacteria
    ## Taxa_00737         28           1228 d__Bacteria
    ## Taxa_00738         25            487 d__Bacteria
    ## Taxa_00739         21            452 d__Bacteria
    ## Taxa_00740         13            233 d__Bacteria
    ## Taxa_00741         27            694 d__Bacteria
    ## Taxa_00742         20            342 d__Bacteria
    ## Taxa_00743         20            718 d__Bacteria
    ## Taxa_00744         12             17 d__Bacteria
    ## Taxa_00745         13             41 d__Bacteria
    ## Taxa_00746          6             10 d__Bacteria
    ## Taxa_00747         27            861 d__Bacteria
    ## Taxa_00748         15             40 d__Bacteria
    ## Taxa_00749          2              2 d__Bacteria
    ## Taxa_00750          7             10 d__Bacteria
    ## Taxa_00751         28            307 d__Bacteria
    ## Taxa_00752         23            142 d__Bacteria
    ## Taxa_00753          1              2 d__Bacteria
    ## Taxa_00754         26            772 d__Bacteria
    ## Taxa_00755          2              6 d__Bacteria
    ## Taxa_00756         13            256 d__Bacteria
    ## Taxa_00757         22            123 d__Bacteria
    ## Taxa_00758          1              1 d__Bacteria
    ## Taxa_00759          1              1 d__Bacteria
    ## Taxa_00760          4              4 d__Bacteria
    ## Taxa_00761          1              1 d__Bacteria
    ## Taxa_00762          1              1 d__Bacteria
    ## Taxa_00763          1              1 d__Bacteria
    ## Taxa_00764          1              1 d__Bacteria
    ## Taxa_00765         18            121 d__Bacteria
    ## Taxa_00766          1              1 d__Bacteria
    ## Taxa_00767         28           5897 d__Bacteria
    ## Taxa_00768          7              7 d__Bacteria
    ## Taxa_00769         15             45 d__Bacteria
    ## Taxa_00770         22            197 d__Bacteria
    ## Taxa_00771          1              1 d__Bacteria
    ## Taxa_00772         12             69 d__Bacteria
    ## Taxa_00773          1              1 d__Bacteria
    ## Taxa_00774         27            368 d__Bacteria
    ## Taxa_00775         21             59 d__Bacteria
    ## Taxa_00776          2              4 d__Bacteria
    ## Taxa_00777          2              2 d__Bacteria
    ## Taxa_00778          3              6 d__Bacteria
    ## Taxa_00779         28           7464 d__Bacteria
    ## Taxa_00780         23             50 d__Bacteria
    ## Taxa_00781          9             11 d__Bacteria
    ## Taxa_00782         15            100 d__Bacteria
    ## Taxa_00783         12             98 d__Bacteria
    ## Taxa_00784          1              1 d__Bacteria
    ## Taxa_00785         26            771 d__Bacteria
    ## Taxa_00786         24            105 d__Bacteria
    ## Taxa_00787          1              1 d__Bacteria
    ## Taxa_00788         15             41 d__Bacteria
    ## Taxa_00789          8             12 d__Bacteria
    ## Taxa_00790         16             77 d__Bacteria
    ## Taxa_00791         10             15 d__Bacteria
    ## Taxa_00792         22             67 d__Bacteria
    ## Taxa_00793          1              1 d__Bacteria
    ## Taxa_00794          2              2 d__Bacteria
    ## Taxa_00795          1              1 d__Bacteria
    ## Taxa_00796         10             13 d__Bacteria
    ## Taxa_00797          7             15 d__Bacteria
    ## Taxa_00798          8             84 d__Bacteria
    ## Taxa_00799          1              1 d__Bacteria
    ## Taxa_00800         28           3234 d__Bacteria
    ## Taxa_00801         28            528 d__Bacteria
    ## Taxa_00802         28            902 d__Bacteria
    ## Taxa_00803         13             27 d__Bacteria
    ## Taxa_00804          1              1 d__Bacteria
    ## Taxa_00805         25            146 d__Bacteria
    ## Taxa_00806          1              1 d__Bacteria
    ## Taxa_00807          8             12 d__Bacteria
    ## Taxa_00808          3             13 d__Bacteria
    ## Taxa_00809         28            714 d__Bacteria
    ## Taxa_00810         28           1143 d__Bacteria
    ## Taxa_00811         16             30 d__Bacteria
    ## Taxa_00812         12             24 d__Bacteria
    ## Taxa_00813          1              1 d__Bacteria
    ## Taxa_00814         13             42 d__Bacteria
    ## Taxa_00815          8             12 d__Bacteria
    ## Taxa_00816          1              1 d__Bacteria
    ## Taxa_00817          7              8 d__Bacteria
    ## Taxa_00818         14             24 d__Bacteria
    ## Taxa_00819         28            635 d__Bacteria
    ## Taxa_00820         27            325 d__Bacteria
    ## Taxa_00821         25            144 d__Bacteria
    ## Taxa_00822         26            152 d__Bacteria
    ## Taxa_00823          3              3 d__Bacteria
    ## Taxa_00824         25            137 d__Bacteria
    ## Taxa_00825         23             65 d__Bacteria
    ## Taxa_00826         15             30 d__Bacteria
    ## Taxa_00827          3              4 d__Bacteria
    ## Taxa_00828          3             25 d__Bacteria
    ## Taxa_00829          9             19 d__Bacteria
    ## Taxa_00830          3              4 d__Bacteria
    ## Taxa_00831          2              2 d__Bacteria
    ## Taxa_00832         28           8828 d__Bacteria
    ## Taxa_00833         27            465 d__Bacteria
    ## Taxa_00834         19             53 d__Bacteria
    ## Taxa_00835          2              2 d__Bacteria
    ## Taxa_00836         14             19 d__Bacteria
    ## Taxa_00837          1              1 d__Bacteria
    ## Taxa_00838          1              1 d__Bacteria
    ## Taxa_00839          3              4 d__Bacteria
    ## Taxa_00840         28            318 d__Bacteria
    ## Taxa_00841          8             35 d__Bacteria
    ## Taxa_00842         16             25 d__Bacteria
    ## Taxa_00843          3              5 d__Bacteria
    ## Taxa_00844          1              1 d__Bacteria
    ## Taxa_00845          3              6 d__Bacteria
    ## Taxa_00846          2              9 d__Bacteria
    ## Taxa_00847          1              1 d__Bacteria
    ## Taxa_00848          6              6 d__Bacteria
    ## Taxa_00849          2             12 d__Bacteria
    ## Taxa_00850          1              1 d__Bacteria
    ## Taxa_00851          1              1 d__Bacteria
    ## Taxa_00852          3             13 d__Bacteria
    ## Taxa_00853          8             13 d__Bacteria
    ## Taxa_00854         28            503 d__Bacteria
    ## Taxa_00855          3              9 d__Bacteria
    ## Taxa_00856          1              1 d__Bacteria
    ## Taxa_00857          3              3 d__Bacteria
    ## Taxa_00858         22             78 d__Bacteria
    ## Taxa_00859         28            980 d__Bacteria
    ## Taxa_00860         27            104 d__Bacteria
    ## Taxa_00861         27            245 d__Bacteria
    ## Taxa_00862         28            615 d__Bacteria
    ## Taxa_00863         18             38 d__Bacteria
    ## Taxa_00864         13             18 d__Bacteria
    ## Taxa_00865          8             14 d__Bacteria
    ## Taxa_00866          1              1 d__Bacteria
    ## Taxa_00867          2              3 d__Bacteria
    ## Taxa_00868         14             22 d__Bacteria
    ## Taxa_00869          1              1 d__Bacteria
    ## Taxa_00870          1              1 d__Bacteria
    ## Taxa_00871          8             21 d__Bacteria
    ## Taxa_00872          1              1 d__Bacteria
    ## Taxa_00873          1              3 d__Bacteria
    ## Taxa_00874          1              1 d__Bacteria
    ## Taxa_00875         18            207 d__Bacteria
    ## Taxa_00876          1              1 d__Bacteria
    ## Taxa_00877         10             23 d__Bacteria
    ## Taxa_00878         21            645 d__Bacteria
    ## Taxa_00879         28           1145 d__Bacteria
    ## Taxa_00880         13            235 d__Bacteria
    ## Taxa_00881          4              8 d__Bacteria
    ## Taxa_00882         20            222 d__Bacteria
    ## Taxa_00883         26            177 d__Bacteria
    ## Taxa_00884          1              1 d__Bacteria
    ## Taxa_00885         28           3290 d__Bacteria
    ## Taxa_00886         17             80 d__Bacteria
    ## Taxa_00887         10             35 d__Bacteria
    ## Taxa_00888         16             52 d__Bacteria
    ## Taxa_00889         19            198 d__Bacteria
    ## Taxa_00890         28          10779 d__Bacteria
    ## Taxa_00891          8             77 d__Bacteria
    ## Taxa_00892         27            278 d__Bacteria
    ## Taxa_00893         28           1515 d__Bacteria
    ## Taxa_00894          2              2 d__Bacteria
    ## Taxa_00895         24            363 d__Bacteria
    ## Taxa_00896         17            117 d__Bacteria
    ## Taxa_00897         15             97 d__Bacteria
    ## Taxa_00898         17            143 d__Bacteria
    ## Taxa_00899          1              2 d__Bacteria
    ## Taxa_00900          1              2 d__Bacteria
    ## Taxa_00901          6             24 d__Bacteria
    ## Taxa_00902          2              3 d__Bacteria
    ## Taxa_00903         15             67 d__Bacteria
    ## Taxa_00904          8             17 d__Bacteria
    ## Taxa_00905          2            331 d__Bacteria
    ## Taxa_00906          1              1 d__Bacteria
    ## Taxa_00907          3             18 d__Bacteria
    ## Taxa_00908         27            984 d__Bacteria
    ## Taxa_00909          2             21 d__Bacteria
    ## Taxa_00910          2              2 d__Bacteria
    ## Taxa_00911          8             22 d__Bacteria
    ## Taxa_00912          1              1 d__Bacteria
    ## Taxa_00913          4              9 d__Bacteria
    ## Taxa_00914         14             63 d__Bacteria
    ## Taxa_00915          5              5 d__Bacteria
    ## Taxa_00916          3              4 d__Bacteria
    ## Taxa_00917          6             14 d__Bacteria
    ## Taxa_00918          7              9 d__Bacteria
    ## Taxa_00919          6             16 d__Bacteria
    ## Taxa_00920          2              3 d__Bacteria
    ## Taxa_00921          7             12 d__Bacteria
    ## Taxa_00922          5             22 d__Bacteria
    ## Taxa_00923          7             15 d__Bacteria
    ## Taxa_00924          1              1 d__Bacteria
    ## Taxa_00925          1              1 d__Bacteria
    ## Taxa_00926          1              1 d__Bacteria
    ## Taxa_00927         28           4614 d__Bacteria
    ## Taxa_00928         15             20 d__Bacteria
    ## Taxa_00929         24             96 d__Bacteria
    ## Taxa_00930          5              6 d__Bacteria
    ## Taxa_00931         28           4499 d__Bacteria
    ## Taxa_00932         10             21 d__Bacteria
    ## Taxa_00933          1              1 d__Bacteria
    ## Taxa_00934         15             53 d__Bacteria
    ## Taxa_00935         28           2359 d__Bacteria
    ## Taxa_00936         28          25850 d__Bacteria
    ## Taxa_00937         28           3677 d__Bacteria
    ## Taxa_00938         27            313 d__Bacteria
    ## Taxa_00939          9             10 d__Bacteria
    ## Taxa_00940         28          11042 d__Bacteria
    ## Taxa_00941         28           1049 d__Bacteria
    ## Taxa_00942         28            274 d__Bacteria
    ## Taxa_00943          3             18 d__Bacteria
    ## Taxa_00944         13            112 d__Bacteria
    ## Taxa_00945         13             18 d__Bacteria
    ## Taxa_00946         17             79 d__Bacteria
    ## Taxa_00947         25            397 d__Bacteria
    ## Taxa_00948         23            240 d__Bacteria
    ## Taxa_00949         28            773  d__unknown
    ##                                    Phylum
    ## Taxa_00000                           <NA>
    ## Taxa_00001               p__Crenarchaeota
    ## Taxa_00002                           <NA>
    ## Taxa_00003               p__Acidobacteria
    ## Taxa_00004               p__Acidobacteria
    ## Taxa_00005               p__Acidobacteria
    ## Taxa_00006               p__Acidobacteria
    ## Taxa_00007               p__Acidobacteria
    ## Taxa_00008               p__Acidobacteria
    ## Taxa_00009               p__Acidobacteria
    ## Taxa_00010               p__Acidobacteria
    ## Taxa_00011               p__Acidobacteria
    ## Taxa_00012               p__Acidobacteria
    ## Taxa_00013               p__Acidobacteria
    ## Taxa_00014               p__Acidobacteria
    ## Taxa_00015               p__Acidobacteria
    ## Taxa_00016               p__Acidobacteria
    ## Taxa_00017               p__Acidobacteria
    ## Taxa_00018               p__Acidobacteria
    ## Taxa_00019               p__Acidobacteria
    ## Taxa_00020               p__Acidobacteria
    ## Taxa_00021               p__Acidobacteria
    ## Taxa_00022               p__Acidobacteria
    ## Taxa_00023               p__Acidobacteria
    ## Taxa_00024               p__Acidobacteria
    ## Taxa_00025               p__Acidobacteria
    ## Taxa_00026               p__Acidobacteria
    ## Taxa_00027               p__Acidobacteria
    ## Taxa_00028               p__Acidobacteria
    ## Taxa_00029               p__Acidobacteria
    ## Taxa_00030               p__Acidobacteria
    ## Taxa_00031               p__Acidobacteria
    ## Taxa_00032               p__Acidobacteria
    ## Taxa_00033               p__Acidobacteria
    ## Taxa_00034               p__Acidobacteria
    ## Taxa_00035               p__Acidobacteria
    ## Taxa_00036               p__Acidobacteria
    ## Taxa_00037               p__Acidobacteria
    ## Taxa_00038               p__Acidobacteria
    ## Taxa_00039               p__Acidobacteria
    ## Taxa_00040               p__Acidobacteria
    ## Taxa_00041               p__Acidobacteria
    ## Taxa_00042               p__Acidobacteria
    ## Taxa_00043               p__Acidobacteria
    ## Taxa_00044              p__Actinobacteria
    ## Taxa_00045              p__Actinobacteria
    ## Taxa_00046              p__Actinobacteria
    ## Taxa_00047              p__Actinobacteria
    ## Taxa_00048              p__Actinobacteria
    ## Taxa_00049              p__Actinobacteria
    ## Taxa_00050              p__Actinobacteria
    ## Taxa_00051              p__Actinobacteria
    ## Taxa_00052              p__Actinobacteria
    ## Taxa_00053              p__Actinobacteria
    ## Taxa_00054              p__Actinobacteria
    ## Taxa_00055              p__Actinobacteria
    ## Taxa_00056              p__Actinobacteria
    ## Taxa_00057              p__Actinobacteria
    ## Taxa_00058              p__Actinobacteria
    ## Taxa_00059              p__Actinobacteria
    ## Taxa_00060              p__Actinobacteria
    ## Taxa_00061              p__Actinobacteria
    ## Taxa_00062              p__Actinobacteria
    ## Taxa_00063              p__Actinobacteria
    ## Taxa_00064              p__Actinobacteria
    ## Taxa_00065              p__Actinobacteria
    ## Taxa_00066              p__Actinobacteria
    ## Taxa_00067              p__Actinobacteria
    ## Taxa_00068              p__Actinobacteria
    ## Taxa_00069              p__Actinobacteria
    ## Taxa_00070              p__Actinobacteria
    ## Taxa_00071              p__Actinobacteria
    ## Taxa_00072              p__Actinobacteria
    ## Taxa_00073              p__Actinobacteria
    ## Taxa_00074              p__Actinobacteria
    ## Taxa_00075              p__Actinobacteria
    ## Taxa_00076              p__Actinobacteria
    ## Taxa_00077              p__Actinobacteria
    ## Taxa_00078              p__Actinobacteria
    ## Taxa_00079              p__Actinobacteria
    ## Taxa_00080              p__Actinobacteria
    ## Taxa_00081              p__Actinobacteria
    ## Taxa_00082              p__Actinobacteria
    ## Taxa_00083              p__Actinobacteria
    ## Taxa_00084              p__Actinobacteria
    ## Taxa_00085              p__Actinobacteria
    ## Taxa_00086              p__Actinobacteria
    ## Taxa_00087              p__Actinobacteria
    ## Taxa_00088              p__Actinobacteria
    ## Taxa_00089              p__Actinobacteria
    ## Taxa_00090              p__Actinobacteria
    ## Taxa_00091              p__Actinobacteria
    ## Taxa_00092              p__Actinobacteria
    ## Taxa_00093              p__Actinobacteria
    ## Taxa_00094              p__Actinobacteria
    ## Taxa_00095              p__Actinobacteria
    ## Taxa_00096              p__Actinobacteria
    ## Taxa_00097              p__Actinobacteria
    ## Taxa_00098              p__Actinobacteria
    ## Taxa_00099              p__Actinobacteria
    ## Taxa_00100              p__Actinobacteria
    ## Taxa_00101              p__Actinobacteria
    ## Taxa_00102              p__Actinobacteria
    ## Taxa_00103              p__Actinobacteria
    ## Taxa_00104              p__Actinobacteria
    ## Taxa_00105              p__Actinobacteria
    ## Taxa_00106              p__Actinobacteria
    ## Taxa_00107              p__Actinobacteria
    ## Taxa_00108              p__Actinobacteria
    ## Taxa_00109              p__Actinobacteria
    ## Taxa_00110              p__Actinobacteria
    ## Taxa_00111              p__Actinobacteria
    ## Taxa_00112              p__Actinobacteria
    ## Taxa_00113              p__Actinobacteria
    ## Taxa_00114              p__Actinobacteria
    ## Taxa_00115              p__Actinobacteria
    ## Taxa_00116              p__Actinobacteria
    ## Taxa_00117              p__Actinobacteria
    ## Taxa_00118              p__Actinobacteria
    ## Taxa_00119              p__Actinobacteria
    ## Taxa_00120              p__Actinobacteria
    ## Taxa_00121              p__Actinobacteria
    ## Taxa_00122              p__Actinobacteria
    ## Taxa_00123              p__Actinobacteria
    ## Taxa_00124              p__Actinobacteria
    ## Taxa_00125              p__Actinobacteria
    ## Taxa_00126              p__Actinobacteria
    ## Taxa_00127              p__Actinobacteria
    ## Taxa_00128              p__Actinobacteria
    ## Taxa_00129              p__Actinobacteria
    ## Taxa_00130              p__Actinobacteria
    ## Taxa_00131              p__Actinobacteria
    ## Taxa_00132              p__Actinobacteria
    ## Taxa_00133              p__Actinobacteria
    ## Taxa_00134              p__Actinobacteria
    ## Taxa_00135              p__Actinobacteria
    ## Taxa_00136              p__Actinobacteria
    ## Taxa_00137              p__Actinobacteria
    ## Taxa_00138              p__Actinobacteria
    ## Taxa_00139              p__Actinobacteria
    ## Taxa_00140              p__Actinobacteria
    ## Taxa_00141              p__Actinobacteria
    ## Taxa_00142              p__Actinobacteria
    ## Taxa_00143              p__Actinobacteria
    ## Taxa_00144              p__Actinobacteria
    ## Taxa_00145              p__Actinobacteria
    ## Taxa_00146              p__Actinobacteria
    ## Taxa_00147              p__Actinobacteria
    ## Taxa_00148              p__Actinobacteria
    ## Taxa_00149              p__Actinobacteria
    ## Taxa_00150              p__Actinobacteria
    ## Taxa_00151              p__Actinobacteria
    ## Taxa_00152              p__Actinobacteria
    ## Taxa_00153              p__Actinobacteria
    ## Taxa_00154              p__Actinobacteria
    ## Taxa_00155              p__Actinobacteria
    ## Taxa_00156              p__Actinobacteria
    ## Taxa_00157              p__Actinobacteria
    ## Taxa_00158              p__Actinobacteria
    ## Taxa_00159              p__Actinobacteria
    ## Taxa_00160              p__Actinobacteria
    ## Taxa_00161              p__Actinobacteria
    ## Taxa_00162              p__Actinobacteria
    ## Taxa_00163              p__Actinobacteria
    ## Taxa_00164              p__Actinobacteria
    ## Taxa_00165              p__Actinobacteria
    ## Taxa_00166              p__Actinobacteria
    ## Taxa_00167              p__Actinobacteria
    ## Taxa_00168              p__Actinobacteria
    ## Taxa_00169              p__Actinobacteria
    ## Taxa_00170              p__Actinobacteria
    ## Taxa_00171              p__Actinobacteria
    ## Taxa_00172              p__Actinobacteria
    ## Taxa_00173              p__Actinobacteria
    ## Taxa_00174              p__Actinobacteria
    ## Taxa_00175              p__Actinobacteria
    ## Taxa_00176              p__Actinobacteria
    ## Taxa_00177              p__Actinobacteria
    ## Taxa_00178              p__Actinobacteria
    ## Taxa_00179              p__Actinobacteria
    ## Taxa_00180              p__Actinobacteria
    ## Taxa_00181              p__Actinobacteria
    ## Taxa_00182              p__Actinobacteria
    ## Taxa_00183              p__Actinobacteria
    ## Taxa_00184              p__Actinobacteria
    ## Taxa_00185              p__Actinobacteria
    ## Taxa_00186              p__Actinobacteria
    ## Taxa_00187              p__Actinobacteria
    ## Taxa_00188              p__Actinobacteria
    ## Taxa_00189              p__Actinobacteria
    ## Taxa_00190              p__Actinobacteria
    ## Taxa_00191              p__Actinobacteria
    ## Taxa_00192              p__Actinobacteria
    ## Taxa_00193              p__Actinobacteria
    ## Taxa_00194              p__Actinobacteria
    ## Taxa_00195              p__Actinobacteria
    ## Taxa_00196              p__Actinobacteria
    ## Taxa_00197              p__Actinobacteria
    ## Taxa_00198              p__Actinobacteria
    ## Taxa_00199              p__Actinobacteria
    ## Taxa_00200              p__Actinobacteria
    ## Taxa_00201              p__Actinobacteria
    ## Taxa_00202              p__Actinobacteria
    ## Taxa_00203              p__Actinobacteria
    ## Taxa_00204              p__Actinobacteria
    ## Taxa_00205              p__Actinobacteria
    ## Taxa_00206              p__Actinobacteria
    ## Taxa_00207              p__Actinobacteria
    ## Taxa_00208              p__Actinobacteria
    ## Taxa_00209              p__Actinobacteria
    ## Taxa_00210              p__Actinobacteria
    ## Taxa_00211              p__Actinobacteria
    ## Taxa_00212              p__Actinobacteria
    ## Taxa_00213              p__Actinobacteria
    ## Taxa_00214              p__Actinobacteria
    ## Taxa_00215              p__Actinobacteria
    ## Taxa_00216              p__Actinobacteria
    ## Taxa_00217              p__Actinobacteria
    ## Taxa_00218                   p__Aquificae
    ## Taxa_00219                   p__Aquificae
    ## Taxa_00220                   p__Aquificae
    ## Taxa_00221             p__Armatimonadetes
    ## Taxa_00222             p__Armatimonadetes
    ## Taxa_00223             p__Armatimonadetes
    ## Taxa_00224             p__Armatimonadetes
    ## Taxa_00225             p__Armatimonadetes
    ## Taxa_00226             p__Armatimonadetes
    ## Taxa_00227             p__Armatimonadetes
    ## Taxa_00228               p__Bacteroidetes
    ## Taxa_00229               p__Bacteroidetes
    ## Taxa_00230               p__Bacteroidetes
    ## Taxa_00231               p__Bacteroidetes
    ## Taxa_00232               p__Bacteroidetes
    ## Taxa_00233               p__Bacteroidetes
    ## Taxa_00234               p__Bacteroidetes
    ## Taxa_00235               p__Bacteroidetes
    ## Taxa_00236               p__Bacteroidetes
    ## Taxa_00237               p__Bacteroidetes
    ## Taxa_00238               p__Bacteroidetes
    ## Taxa_00239               p__Bacteroidetes
    ## Taxa_00240               p__Bacteroidetes
    ## Taxa_00241               p__Bacteroidetes
    ## Taxa_00242               p__Bacteroidetes
    ## Taxa_00243               p__Bacteroidetes
    ## Taxa_00244               p__Bacteroidetes
    ## Taxa_00245               p__Bacteroidetes
    ## Taxa_00246               p__Bacteroidetes
    ## Taxa_00247               p__Bacteroidetes
    ## Taxa_00248               p__Bacteroidetes
    ## Taxa_00249               p__Bacteroidetes
    ## Taxa_00250               p__Bacteroidetes
    ## Taxa_00251               p__Bacteroidetes
    ## Taxa_00252               p__Bacteroidetes
    ## Taxa_00253               p__Bacteroidetes
    ## Taxa_00254               p__Bacteroidetes
    ## Taxa_00255               p__Bacteroidetes
    ## Taxa_00256               p__Bacteroidetes
    ## Taxa_00257               p__Bacteroidetes
    ## Taxa_00258               p__Bacteroidetes
    ## Taxa_00259               p__Bacteroidetes
    ## Taxa_00260               p__Bacteroidetes
    ## Taxa_00261               p__Bacteroidetes
    ## Taxa_00262               p__Bacteroidetes
    ## Taxa_00263               p__Bacteroidetes
    ## Taxa_00264               p__Bacteroidetes
    ## Taxa_00265               p__Bacteroidetes
    ## Taxa_00266               p__Bacteroidetes
    ## Taxa_00267               p__Bacteroidetes
    ## Taxa_00268               p__Bacteroidetes
    ## Taxa_00269               p__Bacteroidetes
    ## Taxa_00270               p__Bacteroidetes
    ## Taxa_00271               p__Bacteroidetes
    ## Taxa_00272               p__Bacteroidetes
    ## Taxa_00273               p__Bacteroidetes
    ## Taxa_00274               p__Bacteroidetes
    ## Taxa_00275               p__Bacteroidetes
    ## Taxa_00276               p__Bacteroidetes
    ## Taxa_00277               p__Bacteroidetes
    ## Taxa_00278               p__Bacteroidetes
    ## Taxa_00279               p__Bacteroidetes
    ## Taxa_00280               p__Bacteroidetes
    ## Taxa_00281               p__Bacteroidetes
    ## Taxa_00282               p__Bacteroidetes
    ## Taxa_00283               p__Bacteroidetes
    ## Taxa_00284               p__Bacteroidetes
    ## Taxa_00285               p__Bacteroidetes
    ## Taxa_00286               p__Bacteroidetes
    ## Taxa_00287               p__Bacteroidetes
    ## Taxa_00288               p__Bacteroidetes
    ## Taxa_00289               p__Bacteroidetes
    ## Taxa_00290               p__Bacteroidetes
    ## Taxa_00291               p__Bacteroidetes
    ## Taxa_00292               p__Bacteroidetes
    ## Taxa_00293               p__Bacteroidetes
    ## Taxa_00294               p__Bacteroidetes
    ## Taxa_00295               p__Bacteroidetes
    ## Taxa_00296               p__Bacteroidetes
    ## Taxa_00297               p__Bacteroidetes
    ## Taxa_00298               p__Bacteroidetes
    ## Taxa_00299               p__Bacteroidetes
    ## Taxa_00300               p__Bacteroidetes
    ## Taxa_00301               p__Bacteroidetes
    ## Taxa_00302               p__Bacteroidetes
    ## Taxa_00303               p__Bacteroidetes
    ## Taxa_00304                        p__BRC1
    ## Taxa_00305    p__candidate division WPS-1
    ## Taxa_00306    p__candidate division WPS-2
    ## Taxa_00307      p__candidate division ZB3
    ## Taxa_00308 p__Candidatus Saccharibacteria
    ## Taxa_00309                  p__Chlamydiae
    ## Taxa_00310                  p__Chlamydiae
    ## Taxa_00311                  p__Chlamydiae
    ## Taxa_00312                  p__Chlamydiae
    ## Taxa_00313                  p__Chlamydiae
    ## Taxa_00314                 p__Chloroflexi
    ## Taxa_00315                 p__Chloroflexi
    ## Taxa_00316                 p__Chloroflexi
    ## Taxa_00317                 p__Chloroflexi
    ## Taxa_00318                 p__Chloroflexi
    ## Taxa_00319                 p__Chloroflexi
    ## Taxa_00320                 p__Chloroflexi
    ## Taxa_00321                 p__Chloroflexi
    ## Taxa_00322                 p__Chloroflexi
    ## Taxa_00323                 p__Chloroflexi
    ## Taxa_00324                 p__Chloroflexi
    ## Taxa_00325                 p__Chloroflexi
    ## Taxa_00326                 p__Chloroflexi
    ## Taxa_00327                 p__Chloroflexi
    ## Taxa_00328                 p__Chloroflexi
    ## Taxa_00329                 p__Chloroflexi
    ## Taxa_00330                 p__Chloroflexi
    ## Taxa_00331                 p__Chloroflexi
    ## Taxa_00332                 p__Chloroflexi
    ## Taxa_00333                 p__Chloroflexi
    ## Taxa_00334                 p__Chloroflexi
    ## Taxa_00335                 p__Chloroflexi
    ## Taxa_00336                 p__Chloroflexi
    ## Taxa_00337                 p__Chloroflexi
    ## Taxa_00338                 p__Chloroflexi
    ## Taxa_00339                 p__Chloroflexi
    ## Taxa_00340                 p__Chloroflexi
    ## Taxa_00341   p__Cyanobacteria/Chloroplast
    ## Taxa_00342   p__Cyanobacteria/Chloroplast
    ## Taxa_00343   p__Cyanobacteria/Chloroplast
    ## Taxa_00344   p__Cyanobacteria/Chloroplast
    ## Taxa_00345   p__Cyanobacteria/Chloroplast
    ## Taxa_00346   p__Cyanobacteria/Chloroplast
    ## Taxa_00347   p__Cyanobacteria/Chloroplast
    ## Taxa_00348   p__Cyanobacteria/Chloroplast
    ## Taxa_00349         p__Deinococcus-Thermus
    ## Taxa_00350               p__Elusimicrobia
    ## Taxa_00351               p__Elusimicrobia
    ## Taxa_00352               p__Elusimicrobia
    ## Taxa_00353                  p__Firmicutes
    ## Taxa_00354                  p__Firmicutes
    ## Taxa_00355                  p__Firmicutes
    ## Taxa_00356                  p__Firmicutes
    ## Taxa_00357                  p__Firmicutes
    ## Taxa_00358                  p__Firmicutes
    ## Taxa_00359                  p__Firmicutes
    ## Taxa_00360                  p__Firmicutes
    ## Taxa_00361                  p__Firmicutes
    ## Taxa_00362                  p__Firmicutes
    ## Taxa_00363                  p__Firmicutes
    ## Taxa_00364                  p__Firmicutes
    ## Taxa_00365                  p__Firmicutes
    ## Taxa_00366                  p__Firmicutes
    ## Taxa_00367                  p__Firmicutes
    ## Taxa_00368                  p__Firmicutes
    ## Taxa_00369                  p__Firmicutes
    ## Taxa_00370                  p__Firmicutes
    ## Taxa_00371                  p__Firmicutes
    ## Taxa_00372                  p__Firmicutes
    ## Taxa_00373                  p__Firmicutes
    ## Taxa_00374                  p__Firmicutes
    ## Taxa_00375                  p__Firmicutes
    ## Taxa_00376                  p__Firmicutes
    ## Taxa_00377                  p__Firmicutes
    ## Taxa_00378                  p__Firmicutes
    ## Taxa_00379                  p__Firmicutes
    ## Taxa_00380                  p__Firmicutes
    ## Taxa_00381                  p__Firmicutes
    ## Taxa_00382                  p__Firmicutes
    ## Taxa_00383                  p__Firmicutes
    ## Taxa_00384                  p__Firmicutes
    ## Taxa_00385                  p__Firmicutes
    ## Taxa_00386                  p__Firmicutes
    ## Taxa_00387                  p__Firmicutes
    ## Taxa_00388                  p__Firmicutes
    ## Taxa_00389                  p__Firmicutes
    ## Taxa_00390                  p__Firmicutes
    ## Taxa_00391                  p__Firmicutes
    ## Taxa_00392                  p__Firmicutes
    ## Taxa_00393                  p__Firmicutes
    ## Taxa_00394                  p__Firmicutes
    ## Taxa_00395                  p__Firmicutes
    ## Taxa_00396                  p__Firmicutes
    ## Taxa_00397                  p__Firmicutes
    ## Taxa_00398                  p__Firmicutes
    ## Taxa_00399                  p__Firmicutes
    ## Taxa_00400                  p__Firmicutes
    ## Taxa_00401                  p__Firmicutes
    ## Taxa_00402                  p__Firmicutes
    ## Taxa_00403                  p__Firmicutes
    ## Taxa_00404                  p__Firmicutes
    ## Taxa_00405                  p__Firmicutes
    ## Taxa_00406                  p__Firmicutes
    ## Taxa_00407                  p__Firmicutes
    ## Taxa_00408                  p__Firmicutes
    ## Taxa_00409                  p__Firmicutes
    ## Taxa_00410                  p__Firmicutes
    ## Taxa_00411                  p__Firmicutes
    ## Taxa_00412                  p__Firmicutes
    ## Taxa_00413                  p__Firmicutes
    ## Taxa_00414                  p__Firmicutes
    ## Taxa_00415                  p__Firmicutes
    ## Taxa_00416                  p__Firmicutes
    ## Taxa_00417                  p__Firmicutes
    ## Taxa_00418                  p__Firmicutes
    ## Taxa_00419                  p__Firmicutes
    ## Taxa_00420                  p__Firmicutes
    ## Taxa_00421                  p__Firmicutes
    ## Taxa_00422                  p__Firmicutes
    ## Taxa_00423                  p__Firmicutes
    ## Taxa_00424                  p__Firmicutes
    ## Taxa_00425                  p__Firmicutes
    ## Taxa_00426                  p__Firmicutes
    ## Taxa_00427                  p__Firmicutes
    ## Taxa_00428                  p__Firmicutes
    ## Taxa_00429                  p__Firmicutes
    ## Taxa_00430                  p__Firmicutes
    ## Taxa_00431                  p__Firmicutes
    ## Taxa_00432                  p__Firmicutes
    ## Taxa_00433                  p__Firmicutes
    ## Taxa_00434                  p__Firmicutes
    ## Taxa_00435                  p__Firmicutes
    ## Taxa_00436                  p__Firmicutes
    ## Taxa_00437                  p__Firmicutes
    ## Taxa_00438                  p__Firmicutes
    ## Taxa_00439                  p__Firmicutes
    ## Taxa_00440                  p__Firmicutes
    ## Taxa_00441                  p__Firmicutes
    ## Taxa_00442                  p__Firmicutes
    ## Taxa_00443                  p__Firmicutes
    ## Taxa_00444                  p__Firmicutes
    ## Taxa_00445                  p__Firmicutes
    ## Taxa_00446                  p__Firmicutes
    ## Taxa_00447                  p__Firmicutes
    ## Taxa_00448                  p__Firmicutes
    ## Taxa_00449                  p__Firmicutes
    ## Taxa_00450                  p__Firmicutes
    ## Taxa_00451                  p__Firmicutes
    ## Taxa_00452                  p__Firmicutes
    ## Taxa_00453                  p__Firmicutes
    ## Taxa_00454                  p__Firmicutes
    ## Taxa_00455                  p__Firmicutes
    ## Taxa_00456                  p__Firmicutes
    ## Taxa_00457                  p__Firmicutes
    ## Taxa_00458                  p__Firmicutes
    ## Taxa_00459                  p__Firmicutes
    ## Taxa_00460                  p__Firmicutes
    ## Taxa_00461                  p__Firmicutes
    ## Taxa_00462                  p__Firmicutes
    ## Taxa_00463                  p__Firmicutes
    ## Taxa_00464                  p__Firmicutes
    ## Taxa_00465                  p__Firmicutes
    ## Taxa_00466                  p__Firmicutes
    ## Taxa_00467            p__Gemmatimonadetes
    ## Taxa_00468             p__Hydrogenedentes
    ## Taxa_00469             p__Latescibacteria
    ## Taxa_00470               p__Lentisphaerae
    ## Taxa_00471               p__Lentisphaerae
    ## Taxa_00472              p__Microgenomates
    ## Taxa_00473                 p__Nitrospirae
    ## Taxa_00474                p__Omnitrophica
    ## Taxa_00475               p__Parcubacteria
    ## Taxa_00476              p__Planctomycetes
    ## Taxa_00477              p__Planctomycetes
    ## Taxa_00478              p__Planctomycetes
    ## Taxa_00479              p__Planctomycetes
    ## Taxa_00480              p__Planctomycetes
    ## Taxa_00481              p__Planctomycetes
    ## Taxa_00482              p__Planctomycetes
    ## Taxa_00483              p__Planctomycetes
    ## Taxa_00484              p__Planctomycetes
    ## Taxa_00485              p__Planctomycetes
    ## Taxa_00486              p__Planctomycetes
    ## Taxa_00487              p__Planctomycetes
    ## Taxa_00488              p__Planctomycetes
    ## Taxa_00489              p__Planctomycetes
    ## Taxa_00490              p__Planctomycetes
    ## Taxa_00491              p__Planctomycetes
    ## Taxa_00492              p__Planctomycetes
    ## Taxa_00493              p__Planctomycetes
    ## Taxa_00494              p__Planctomycetes
    ## Taxa_00495              p__Planctomycetes
    ## Taxa_00496              p__Planctomycetes
    ## Taxa_00497              p__Planctomycetes
    ## Taxa_00498              p__Proteobacteria
    ## Taxa_00499              p__Proteobacteria
    ## Taxa_00500              p__Proteobacteria
    ## Taxa_00501              p__Proteobacteria
    ## Taxa_00502              p__Proteobacteria
    ## Taxa_00503              p__Proteobacteria
    ## Taxa_00504              p__Proteobacteria
    ## Taxa_00505              p__Proteobacteria
    ## Taxa_00506              p__Proteobacteria
    ## Taxa_00507              p__Proteobacteria
    ## Taxa_00508              p__Proteobacteria
    ## Taxa_00509              p__Proteobacteria
    ## Taxa_00510              p__Proteobacteria
    ## Taxa_00511              p__Proteobacteria
    ## Taxa_00512              p__Proteobacteria
    ## Taxa_00513              p__Proteobacteria
    ## Taxa_00514              p__Proteobacteria
    ## Taxa_00515              p__Proteobacteria
    ## Taxa_00516              p__Proteobacteria
    ## Taxa_00517              p__Proteobacteria
    ## Taxa_00518              p__Proteobacteria
    ## Taxa_00519              p__Proteobacteria
    ## Taxa_00520              p__Proteobacteria
    ## Taxa_00521              p__Proteobacteria
    ## Taxa_00522              p__Proteobacteria
    ## Taxa_00523              p__Proteobacteria
    ## Taxa_00524              p__Proteobacteria
    ## Taxa_00525              p__Proteobacteria
    ## Taxa_00526              p__Proteobacteria
    ## Taxa_00527              p__Proteobacteria
    ## Taxa_00528              p__Proteobacteria
    ## Taxa_00529              p__Proteobacteria
    ## Taxa_00530              p__Proteobacteria
    ## Taxa_00531              p__Proteobacteria
    ## Taxa_00532              p__Proteobacteria
    ## Taxa_00533              p__Proteobacteria
    ## Taxa_00534              p__Proteobacteria
    ## Taxa_00535              p__Proteobacteria
    ## Taxa_00536              p__Proteobacteria
    ## Taxa_00537              p__Proteobacteria
    ## Taxa_00538              p__Proteobacteria
    ## Taxa_00539              p__Proteobacteria
    ## Taxa_00540              p__Proteobacteria
    ## Taxa_00541              p__Proteobacteria
    ## Taxa_00542              p__Proteobacteria
    ## Taxa_00543              p__Proteobacteria
    ## Taxa_00544              p__Proteobacteria
    ## Taxa_00545              p__Proteobacteria
    ## Taxa_00546              p__Proteobacteria
    ## Taxa_00547              p__Proteobacteria
    ## Taxa_00548              p__Proteobacteria
    ## Taxa_00549              p__Proteobacteria
    ## Taxa_00550              p__Proteobacteria
    ## Taxa_00551              p__Proteobacteria
    ## Taxa_00552              p__Proteobacteria
    ## Taxa_00553              p__Proteobacteria
    ## Taxa_00554              p__Proteobacteria
    ## Taxa_00555              p__Proteobacteria
    ## Taxa_00556              p__Proteobacteria
    ## Taxa_00557              p__Proteobacteria
    ## Taxa_00558              p__Proteobacteria
    ## Taxa_00559              p__Proteobacteria
    ## Taxa_00560              p__Proteobacteria
    ## Taxa_00561              p__Proteobacteria
    ## Taxa_00562              p__Proteobacteria
    ## Taxa_00563              p__Proteobacteria
    ## Taxa_00564              p__Proteobacteria
    ## Taxa_00565              p__Proteobacteria
    ## Taxa_00566              p__Proteobacteria
    ## Taxa_00567              p__Proteobacteria
    ## Taxa_00568              p__Proteobacteria
    ## Taxa_00569              p__Proteobacteria
    ## Taxa_00570              p__Proteobacteria
    ## Taxa_00571              p__Proteobacteria
    ## Taxa_00572              p__Proteobacteria
    ## Taxa_00573              p__Proteobacteria
    ## Taxa_00574              p__Proteobacteria
    ## Taxa_00575              p__Proteobacteria
    ## Taxa_00576              p__Proteobacteria
    ## Taxa_00577              p__Proteobacteria
    ## Taxa_00578              p__Proteobacteria
    ## Taxa_00579              p__Proteobacteria
    ## Taxa_00580              p__Proteobacteria
    ## Taxa_00581              p__Proteobacteria
    ## Taxa_00582              p__Proteobacteria
    ## Taxa_00583              p__Proteobacteria
    ## Taxa_00584              p__Proteobacteria
    ## Taxa_00585              p__Proteobacteria
    ## Taxa_00586              p__Proteobacteria
    ## Taxa_00587              p__Proteobacteria
    ## Taxa_00588              p__Proteobacteria
    ## Taxa_00589              p__Proteobacteria
    ## Taxa_00590              p__Proteobacteria
    ## Taxa_00591              p__Proteobacteria
    ## Taxa_00592              p__Proteobacteria
    ## Taxa_00593              p__Proteobacteria
    ## Taxa_00594              p__Proteobacteria
    ## Taxa_00595              p__Proteobacteria
    ## Taxa_00596              p__Proteobacteria
    ## Taxa_00597              p__Proteobacteria
    ## Taxa_00598              p__Proteobacteria
    ## Taxa_00599              p__Proteobacteria
    ## Taxa_00600              p__Proteobacteria
    ## Taxa_00601              p__Proteobacteria
    ## Taxa_00602              p__Proteobacteria
    ## Taxa_00603              p__Proteobacteria
    ## Taxa_00604              p__Proteobacteria
    ## Taxa_00605              p__Proteobacteria
    ## Taxa_00606              p__Proteobacteria
    ## Taxa_00607              p__Proteobacteria
    ## Taxa_00608              p__Proteobacteria
    ## Taxa_00609              p__Proteobacteria
    ## Taxa_00610              p__Proteobacteria
    ## Taxa_00611              p__Proteobacteria
    ## Taxa_00612              p__Proteobacteria
    ## Taxa_00613              p__Proteobacteria
    ## Taxa_00614              p__Proteobacteria
    ## Taxa_00615              p__Proteobacteria
    ## Taxa_00616              p__Proteobacteria
    ## Taxa_00617              p__Proteobacteria
    ## Taxa_00618              p__Proteobacteria
    ## Taxa_00619              p__Proteobacteria
    ## Taxa_00620              p__Proteobacteria
    ## Taxa_00621              p__Proteobacteria
    ## Taxa_00622              p__Proteobacteria
    ## Taxa_00623              p__Proteobacteria
    ## Taxa_00624              p__Proteobacteria
    ## Taxa_00625              p__Proteobacteria
    ## Taxa_00626              p__Proteobacteria
    ## Taxa_00627              p__Proteobacteria
    ## Taxa_00628              p__Proteobacteria
    ## Taxa_00629              p__Proteobacteria
    ## Taxa_00630              p__Proteobacteria
    ## Taxa_00631              p__Proteobacteria
    ## Taxa_00632              p__Proteobacteria
    ## Taxa_00633              p__Proteobacteria
    ## Taxa_00634              p__Proteobacteria
    ## Taxa_00635              p__Proteobacteria
    ## Taxa_00636              p__Proteobacteria
    ## Taxa_00637              p__Proteobacteria
    ## Taxa_00638              p__Proteobacteria
    ## Taxa_00639              p__Proteobacteria
    ## Taxa_00640              p__Proteobacteria
    ## Taxa_00641              p__Proteobacteria
    ## Taxa_00642              p__Proteobacteria
    ## Taxa_00643              p__Proteobacteria
    ## Taxa_00644              p__Proteobacteria
    ## Taxa_00645              p__Proteobacteria
    ## Taxa_00646              p__Proteobacteria
    ## Taxa_00647              p__Proteobacteria
    ## Taxa_00648              p__Proteobacteria
    ## Taxa_00649              p__Proteobacteria
    ## Taxa_00650              p__Proteobacteria
    ## Taxa_00651              p__Proteobacteria
    ## Taxa_00652              p__Proteobacteria
    ## Taxa_00653              p__Proteobacteria
    ## Taxa_00654              p__Proteobacteria
    ## Taxa_00655              p__Proteobacteria
    ## Taxa_00656              p__Proteobacteria
    ## Taxa_00657              p__Proteobacteria
    ## Taxa_00658              p__Proteobacteria
    ## Taxa_00659              p__Proteobacteria
    ## Taxa_00660              p__Proteobacteria
    ## Taxa_00661              p__Proteobacteria
    ## Taxa_00662              p__Proteobacteria
    ## Taxa_00663              p__Proteobacteria
    ## Taxa_00664              p__Proteobacteria
    ## Taxa_00665              p__Proteobacteria
    ## Taxa_00666              p__Proteobacteria
    ## Taxa_00667              p__Proteobacteria
    ## Taxa_00668              p__Proteobacteria
    ## Taxa_00669              p__Proteobacteria
    ## Taxa_00670              p__Proteobacteria
    ## Taxa_00671              p__Proteobacteria
    ## Taxa_00672              p__Proteobacteria
    ## Taxa_00673              p__Proteobacteria
    ## Taxa_00674              p__Proteobacteria
    ## Taxa_00675              p__Proteobacteria
    ## Taxa_00676              p__Proteobacteria
    ## Taxa_00677              p__Proteobacteria
    ## Taxa_00678              p__Proteobacteria
    ## Taxa_00679              p__Proteobacteria
    ## Taxa_00680              p__Proteobacteria
    ## Taxa_00681              p__Proteobacteria
    ## Taxa_00682              p__Proteobacteria
    ## Taxa_00683              p__Proteobacteria
    ## Taxa_00684              p__Proteobacteria
    ## Taxa_00685              p__Proteobacteria
    ## Taxa_00686              p__Proteobacteria
    ## Taxa_00687              p__Proteobacteria
    ## Taxa_00688              p__Proteobacteria
    ## Taxa_00689              p__Proteobacteria
    ## Taxa_00690              p__Proteobacteria
    ## Taxa_00691              p__Proteobacteria
    ## Taxa_00692              p__Proteobacteria
    ## Taxa_00693              p__Proteobacteria
    ## Taxa_00694              p__Proteobacteria
    ## Taxa_00695              p__Proteobacteria
    ## Taxa_00696              p__Proteobacteria
    ## Taxa_00697              p__Proteobacteria
    ## Taxa_00698              p__Proteobacteria
    ## Taxa_00699              p__Proteobacteria
    ## Taxa_00700              p__Proteobacteria
    ## Taxa_00701              p__Proteobacteria
    ## Taxa_00702              p__Proteobacteria
    ## Taxa_00703              p__Proteobacteria
    ## Taxa_00704              p__Proteobacteria
    ## Taxa_00705              p__Proteobacteria
    ## Taxa_00706              p__Proteobacteria
    ## Taxa_00707              p__Proteobacteria
    ## Taxa_00708              p__Proteobacteria
    ## Taxa_00709              p__Proteobacteria
    ## Taxa_00710              p__Proteobacteria
    ## Taxa_00711              p__Proteobacteria
    ## Taxa_00712              p__Proteobacteria
    ## Taxa_00713              p__Proteobacteria
    ## Taxa_00714              p__Proteobacteria
    ## Taxa_00715              p__Proteobacteria
    ## Taxa_00716              p__Proteobacteria
    ## Taxa_00717              p__Proteobacteria
    ## Taxa_00718              p__Proteobacteria
    ## Taxa_00719              p__Proteobacteria
    ## Taxa_00720              p__Proteobacteria
    ## Taxa_00721              p__Proteobacteria
    ## Taxa_00722              p__Proteobacteria
    ## Taxa_00723              p__Proteobacteria
    ## Taxa_00724              p__Proteobacteria
    ## Taxa_00725              p__Proteobacteria
    ## Taxa_00726              p__Proteobacteria
    ## Taxa_00727              p__Proteobacteria
    ## Taxa_00728              p__Proteobacteria
    ## Taxa_00729              p__Proteobacteria
    ## Taxa_00730              p__Proteobacteria
    ## Taxa_00731              p__Proteobacteria
    ## Taxa_00732              p__Proteobacteria
    ## Taxa_00733              p__Proteobacteria
    ## Taxa_00734              p__Proteobacteria
    ## Taxa_00735              p__Proteobacteria
    ## Taxa_00736              p__Proteobacteria
    ## Taxa_00737              p__Proteobacteria
    ## Taxa_00738              p__Proteobacteria
    ## Taxa_00739              p__Proteobacteria
    ## Taxa_00740              p__Proteobacteria
    ## Taxa_00741              p__Proteobacteria
    ## Taxa_00742              p__Proteobacteria
    ## Taxa_00743              p__Proteobacteria
    ## Taxa_00744              p__Proteobacteria
    ## Taxa_00745              p__Proteobacteria
    ## Taxa_00746              p__Proteobacteria
    ## Taxa_00747              p__Proteobacteria
    ## Taxa_00748              p__Proteobacteria
    ## Taxa_00749              p__Proteobacteria
    ## Taxa_00750              p__Proteobacteria
    ## Taxa_00751              p__Proteobacteria
    ## Taxa_00752              p__Proteobacteria
    ## Taxa_00753              p__Proteobacteria
    ## Taxa_00754              p__Proteobacteria
    ## Taxa_00755              p__Proteobacteria
    ## Taxa_00756              p__Proteobacteria
    ## Taxa_00757              p__Proteobacteria
    ## Taxa_00758              p__Proteobacteria
    ## Taxa_00759              p__Proteobacteria
    ## Taxa_00760              p__Proteobacteria
    ## Taxa_00761              p__Proteobacteria
    ## Taxa_00762              p__Proteobacteria
    ## Taxa_00763              p__Proteobacteria
    ## Taxa_00764              p__Proteobacteria
    ## Taxa_00765              p__Proteobacteria
    ## Taxa_00766              p__Proteobacteria
    ## Taxa_00767              p__Proteobacteria
    ## Taxa_00768              p__Proteobacteria
    ## Taxa_00769              p__Proteobacteria
    ## Taxa_00770              p__Proteobacteria
    ## Taxa_00771              p__Proteobacteria
    ## Taxa_00772              p__Proteobacteria
    ## Taxa_00773              p__Proteobacteria
    ## Taxa_00774              p__Proteobacteria
    ## Taxa_00775              p__Proteobacteria
    ## Taxa_00776              p__Proteobacteria
    ## Taxa_00777              p__Proteobacteria
    ## Taxa_00778              p__Proteobacteria
    ## Taxa_00779              p__Proteobacteria
    ## Taxa_00780              p__Proteobacteria
    ## Taxa_00781              p__Proteobacteria
    ## Taxa_00782              p__Proteobacteria
    ## Taxa_00783              p__Proteobacteria
    ## Taxa_00784              p__Proteobacteria
    ## Taxa_00785              p__Proteobacteria
    ## Taxa_00786              p__Proteobacteria
    ## Taxa_00787              p__Proteobacteria
    ## Taxa_00788              p__Proteobacteria
    ## Taxa_00789              p__Proteobacteria
    ## Taxa_00790              p__Proteobacteria
    ## Taxa_00791              p__Proteobacteria
    ## Taxa_00792              p__Proteobacteria
    ## Taxa_00793              p__Proteobacteria
    ## Taxa_00794              p__Proteobacteria
    ## Taxa_00795              p__Proteobacteria
    ## Taxa_00796              p__Proteobacteria
    ## Taxa_00797              p__Proteobacteria
    ## Taxa_00798              p__Proteobacteria
    ## Taxa_00799              p__Proteobacteria
    ## Taxa_00800              p__Proteobacteria
    ## Taxa_00801              p__Proteobacteria
    ## Taxa_00802              p__Proteobacteria
    ## Taxa_00803              p__Proteobacteria
    ## Taxa_00804              p__Proteobacteria
    ## Taxa_00805              p__Proteobacteria
    ## Taxa_00806              p__Proteobacteria
    ## Taxa_00807              p__Proteobacteria
    ## Taxa_00808              p__Proteobacteria
    ## Taxa_00809              p__Proteobacteria
    ## Taxa_00810              p__Proteobacteria
    ## Taxa_00811              p__Proteobacteria
    ## Taxa_00812              p__Proteobacteria
    ## Taxa_00813              p__Proteobacteria
    ## Taxa_00814              p__Proteobacteria
    ## Taxa_00815              p__Proteobacteria
    ## Taxa_00816              p__Proteobacteria
    ## Taxa_00817              p__Proteobacteria
    ## Taxa_00818              p__Proteobacteria
    ## Taxa_00819              p__Proteobacteria
    ## Taxa_00820              p__Proteobacteria
    ## Taxa_00821              p__Proteobacteria
    ## Taxa_00822              p__Proteobacteria
    ## Taxa_00823              p__Proteobacteria
    ## Taxa_00824              p__Proteobacteria
    ## Taxa_00825              p__Proteobacteria
    ## Taxa_00826              p__Proteobacteria
    ## Taxa_00827              p__Proteobacteria
    ## Taxa_00828              p__Proteobacteria
    ## Taxa_00829              p__Proteobacteria
    ## Taxa_00830              p__Proteobacteria
    ## Taxa_00831              p__Proteobacteria
    ## Taxa_00832              p__Proteobacteria
    ## Taxa_00833              p__Proteobacteria
    ## Taxa_00834              p__Proteobacteria
    ## Taxa_00835              p__Proteobacteria
    ## Taxa_00836              p__Proteobacteria
    ## Taxa_00837              p__Proteobacteria
    ## Taxa_00838              p__Proteobacteria
    ## Taxa_00839              p__Proteobacteria
    ## Taxa_00840              p__Proteobacteria
    ## Taxa_00841              p__Proteobacteria
    ## Taxa_00842              p__Proteobacteria
    ## Taxa_00843              p__Proteobacteria
    ## Taxa_00844              p__Proteobacteria
    ## Taxa_00845              p__Proteobacteria
    ## Taxa_00846              p__Proteobacteria
    ## Taxa_00847              p__Proteobacteria
    ## Taxa_00848              p__Proteobacteria
    ## Taxa_00849              p__Proteobacteria
    ## Taxa_00850              p__Proteobacteria
    ## Taxa_00851              p__Proteobacteria
    ## Taxa_00852              p__Proteobacteria
    ## Taxa_00853              p__Proteobacteria
    ## Taxa_00854              p__Proteobacteria
    ## Taxa_00855              p__Proteobacteria
    ## Taxa_00856              p__Proteobacteria
    ## Taxa_00857              p__Proteobacteria
    ## Taxa_00858              p__Proteobacteria
    ## Taxa_00859              p__Proteobacteria
    ## Taxa_00860              p__Proteobacteria
    ## Taxa_00861              p__Proteobacteria
    ## Taxa_00862              p__Proteobacteria
    ## Taxa_00863              p__Proteobacteria
    ## Taxa_00864              p__Proteobacteria
    ## Taxa_00865              p__Proteobacteria
    ## Taxa_00866              p__Proteobacteria
    ## Taxa_00867              p__Proteobacteria
    ## Taxa_00868              p__Proteobacteria
    ## Taxa_00869              p__Proteobacteria
    ## Taxa_00870              p__Proteobacteria
    ## Taxa_00871              p__Proteobacteria
    ## Taxa_00872              p__Proteobacteria
    ## Taxa_00873              p__Proteobacteria
    ## Taxa_00874              p__Proteobacteria
    ## Taxa_00875              p__Proteobacteria
    ## Taxa_00876              p__Proteobacteria
    ## Taxa_00877              p__Proteobacteria
    ## Taxa_00878              p__Proteobacteria
    ## Taxa_00879              p__Proteobacteria
    ## Taxa_00880              p__Proteobacteria
    ## Taxa_00881              p__Proteobacteria
    ## Taxa_00882              p__Proteobacteria
    ## Taxa_00883              p__Proteobacteria
    ## Taxa_00884              p__Proteobacteria
    ## Taxa_00885              p__Proteobacteria
    ## Taxa_00886              p__Proteobacteria
    ## Taxa_00887              p__Proteobacteria
    ## Taxa_00888              p__Proteobacteria
    ## Taxa_00889              p__Proteobacteria
    ## Taxa_00890              p__Proteobacteria
    ## Taxa_00891              p__Proteobacteria
    ## Taxa_00892              p__Proteobacteria
    ## Taxa_00893              p__Proteobacteria
    ## Taxa_00894              p__Proteobacteria
    ## Taxa_00895              p__Proteobacteria
    ## Taxa_00896              p__Proteobacteria
    ## Taxa_00897              p__Proteobacteria
    ## Taxa_00898              p__Proteobacteria
    ## Taxa_00899              p__Proteobacteria
    ## Taxa_00900              p__Proteobacteria
    ## Taxa_00901              p__Proteobacteria
    ## Taxa_00902              p__Proteobacteria
    ## Taxa_00903              p__Proteobacteria
    ## Taxa_00904              p__Proteobacteria
    ## Taxa_00905              p__Proteobacteria
    ## Taxa_00906              p__Proteobacteria
    ## Taxa_00907              p__Proteobacteria
    ## Taxa_00908              p__Proteobacteria
    ## Taxa_00909              p__Proteobacteria
    ## Taxa_00910              p__Proteobacteria
    ## Taxa_00911              p__Proteobacteria
    ## Taxa_00912              p__Proteobacteria
    ## Taxa_00913              p__Proteobacteria
    ## Taxa_00914              p__Proteobacteria
    ## Taxa_00915                p__Spirochaetes
    ## Taxa_00916                p__Spirochaetes
    ## Taxa_00917                p__Spirochaetes
    ## Taxa_00918                p__Spirochaetes
    ## Taxa_00919                p__Spirochaetes
    ## Taxa_00920                p__Spirochaetes
    ## Taxa_00921                p__Spirochaetes
    ## Taxa_00922                p__Spirochaetes
    ## Taxa_00923                         p__SR1
    ## Taxa_00924                 p__Tenericutes
    ## Taxa_00925                 p__Tenericutes
    ## Taxa_00926       p__Thermodesulfobacteria
    ## Taxa_00927             p__Verrucomicrobia
    ## Taxa_00928             p__Verrucomicrobia
    ## Taxa_00929             p__Verrucomicrobia
    ## Taxa_00930             p__Verrucomicrobia
    ## Taxa_00931             p__Verrucomicrobia
    ## Taxa_00932             p__Verrucomicrobia
    ## Taxa_00933             p__Verrucomicrobia
    ## Taxa_00934             p__Verrucomicrobia
    ## Taxa_00935             p__Verrucomicrobia
    ## Taxa_00936             p__Verrucomicrobia
    ## Taxa_00937             p__Verrucomicrobia
    ## Taxa_00938             p__Verrucomicrobia
    ## Taxa_00939             p__Verrucomicrobia
    ## Taxa_00940             p__Verrucomicrobia
    ## Taxa_00941             p__Verrucomicrobia
    ## Taxa_00942             p__Verrucomicrobia
    ## Taxa_00943             p__Verrucomicrobia
    ## Taxa_00944             p__Verrucomicrobia
    ## Taxa_00945             p__Verrucomicrobia
    ## Taxa_00946             p__Verrucomicrobia
    ## Taxa_00947             p__Verrucomicrobia
    ## Taxa_00948             p__Verrucomicrobia
    ## Taxa_00949                           <NA>
    ##                                                Class
    ## Taxa_00000                                      <NA>
    ## Taxa_00001                           c__Thermoprotei
    ## Taxa_00002                                      <NA>
    ## Taxa_00003                                      <NA>
    ## Taxa_00004                      c__Acidobacteria_Gp1
    ## Taxa_00005                     c__Acidobacteria_Gp10
    ## Taxa_00006                     c__Acidobacteria_Gp11
    ## Taxa_00007                     c__Acidobacteria_Gp12
    ## Taxa_00008                     c__Acidobacteria_Gp13
    ## Taxa_00009                     c__Acidobacteria_Gp15
    ## Taxa_00010                     c__Acidobacteria_Gp16
    ## Taxa_00011                     c__Acidobacteria_Gp17
    ## Taxa_00012                     c__Acidobacteria_Gp18
    ## Taxa_00013                     c__Acidobacteria_Gp19
    ## Taxa_00014                      c__Acidobacteria_Gp1
    ## Taxa_00015                      c__Acidobacteria_Gp1
    ## Taxa_00016                      c__Acidobacteria_Gp1
    ## Taxa_00017                      c__Acidobacteria_Gp1
    ## Taxa_00018                      c__Acidobacteria_Gp1
    ## Taxa_00019                      c__Acidobacteria_Gp1
    ## Taxa_00020                      c__Acidobacteria_Gp1
    ## Taxa_00021                      c__Acidobacteria_Gp1
    ## Taxa_00022                      c__Acidobacteria_Gp1
    ## Taxa_00023                      c__Acidobacteria_Gp1
    ## Taxa_00024                     c__Acidobacteria_Gp20
    ## Taxa_00025                     c__Acidobacteria_Gp22
    ## Taxa_00026                     c__Acidobacteria_Gp25
    ## Taxa_00027                      c__Acidobacteria_Gp2
    ## Taxa_00028                      c__Acidobacteria_Gp3
    ## Taxa_00029                      c__Acidobacteria_Gp3
    ## Taxa_00030                      c__Acidobacteria_Gp3
    ## Taxa_00031                      c__Acidobacteria_Gp3
    ## Taxa_00032                      c__Acidobacteria_Gp3
    ## Taxa_00033                      c__Acidobacteria_Gp4
    ## Taxa_00034                      c__Acidobacteria_Gp4
    ## Taxa_00035                      c__Acidobacteria_Gp4
    ## Taxa_00036                      c__Acidobacteria_Gp4
    ## Taxa_00037                      c__Acidobacteria_Gp4
    ## Taxa_00038                      c__Acidobacteria_Gp5
    ## Taxa_00039                      c__Acidobacteria_Gp6
    ## Taxa_00040                      c__Acidobacteria_Gp7
    ## Taxa_00041                             c__Holophagae
    ## Taxa_00042                             c__Holophagae
    ## Taxa_00043                             c__Holophagae
    ## Taxa_00044                                      <NA>
    ## Taxa_00045                         c__Actinobacteria
    ## Taxa_00046                         c__Actinobacteria
    ## Taxa_00047                         c__Actinobacteria
    ## Taxa_00048                         c__Actinobacteria
    ## Taxa_00049                         c__Actinobacteria
    ## Taxa_00050                         c__Actinobacteria
    ## Taxa_00051                         c__Actinobacteria
    ## Taxa_00052                         c__Actinobacteria
    ## Taxa_00053                         c__Actinobacteria
    ## Taxa_00054                         c__Actinobacteria
    ## Taxa_00055                         c__Actinobacteria
    ## Taxa_00056                         c__Actinobacteria
    ## Taxa_00057                         c__Actinobacteria
    ## Taxa_00058                         c__Actinobacteria
    ## Taxa_00059                         c__Actinobacteria
    ## Taxa_00060                         c__Actinobacteria
    ## Taxa_00061                         c__Actinobacteria
    ## Taxa_00062                         c__Actinobacteria
    ## Taxa_00063                         c__Actinobacteria
    ## Taxa_00064                         c__Actinobacteria
    ## Taxa_00065                         c__Actinobacteria
    ## Taxa_00066                         c__Actinobacteria
    ## Taxa_00067                         c__Actinobacteria
    ## Taxa_00068                         c__Actinobacteria
    ## Taxa_00069                         c__Actinobacteria
    ## Taxa_00070                         c__Actinobacteria
    ## Taxa_00071                         c__Actinobacteria
    ## Taxa_00072                         c__Actinobacteria
    ## Taxa_00073                         c__Actinobacteria
    ## Taxa_00074                         c__Actinobacteria
    ## Taxa_00075                         c__Actinobacteria
    ## Taxa_00076                         c__Actinobacteria
    ## Taxa_00077                         c__Actinobacteria
    ## Taxa_00078                         c__Actinobacteria
    ## Taxa_00079                         c__Actinobacteria
    ## Taxa_00080                         c__Actinobacteria
    ## Taxa_00081                         c__Actinobacteria
    ## Taxa_00082                         c__Actinobacteria
    ## Taxa_00083                         c__Actinobacteria
    ## Taxa_00084                         c__Actinobacteria
    ## Taxa_00085                         c__Actinobacteria
    ## Taxa_00086                         c__Actinobacteria
    ## Taxa_00087                         c__Actinobacteria
    ## Taxa_00088                         c__Actinobacteria
    ## Taxa_00089                         c__Actinobacteria
    ## Taxa_00090                         c__Actinobacteria
    ## Taxa_00091                         c__Actinobacteria
    ## Taxa_00092                         c__Actinobacteria
    ## Taxa_00093                         c__Actinobacteria
    ## Taxa_00094                         c__Actinobacteria
    ## Taxa_00095                         c__Actinobacteria
    ## Taxa_00096                         c__Actinobacteria
    ## Taxa_00097                         c__Actinobacteria
    ## Taxa_00098                         c__Actinobacteria
    ## Taxa_00099                         c__Actinobacteria
    ## Taxa_00100                         c__Actinobacteria
    ## Taxa_00101                         c__Actinobacteria
    ## Taxa_00102                         c__Actinobacteria
    ## Taxa_00103                         c__Actinobacteria
    ## Taxa_00104                         c__Actinobacteria
    ## Taxa_00105                         c__Actinobacteria
    ## Taxa_00106                         c__Actinobacteria
    ## Taxa_00107                         c__Actinobacteria
    ## Taxa_00108                         c__Actinobacteria
    ## Taxa_00109                         c__Actinobacteria
    ## Taxa_00110                         c__Actinobacteria
    ## Taxa_00111                         c__Actinobacteria
    ## Taxa_00112                         c__Actinobacteria
    ## Taxa_00113                         c__Actinobacteria
    ## Taxa_00114                         c__Actinobacteria
    ## Taxa_00115                         c__Actinobacteria
    ## Taxa_00116                         c__Actinobacteria
    ## Taxa_00117                         c__Actinobacteria
    ## Taxa_00118                         c__Actinobacteria
    ## Taxa_00119                         c__Actinobacteria
    ## Taxa_00120                         c__Actinobacteria
    ## Taxa_00121                         c__Actinobacteria
    ## Taxa_00122                         c__Actinobacteria
    ## Taxa_00123                         c__Actinobacteria
    ## Taxa_00124                         c__Actinobacteria
    ## Taxa_00125                         c__Actinobacteria
    ## Taxa_00126                         c__Actinobacteria
    ## Taxa_00127                         c__Actinobacteria
    ## Taxa_00128                         c__Actinobacteria
    ## Taxa_00129                         c__Actinobacteria
    ## Taxa_00130                         c__Actinobacteria
    ## Taxa_00131                         c__Actinobacteria
    ## Taxa_00132                         c__Actinobacteria
    ## Taxa_00133                         c__Actinobacteria
    ## Taxa_00134                         c__Actinobacteria
    ## Taxa_00135                         c__Actinobacteria
    ## Taxa_00136                         c__Actinobacteria
    ## Taxa_00137                         c__Actinobacteria
    ## Taxa_00138                         c__Actinobacteria
    ## Taxa_00139                         c__Actinobacteria
    ## Taxa_00140                         c__Actinobacteria
    ## Taxa_00141                         c__Actinobacteria
    ## Taxa_00142                         c__Actinobacteria
    ## Taxa_00143                         c__Actinobacteria
    ## Taxa_00144                         c__Actinobacteria
    ## Taxa_00145                         c__Actinobacteria
    ## Taxa_00146                         c__Actinobacteria
    ## Taxa_00147                         c__Actinobacteria
    ## Taxa_00148                         c__Actinobacteria
    ## Taxa_00149                         c__Actinobacteria
    ## Taxa_00150                         c__Actinobacteria
    ## Taxa_00151                         c__Actinobacteria
    ## Taxa_00152                         c__Actinobacteria
    ## Taxa_00153                         c__Actinobacteria
    ## Taxa_00154                         c__Actinobacteria
    ## Taxa_00155                         c__Actinobacteria
    ## Taxa_00156                         c__Actinobacteria
    ## Taxa_00157                         c__Actinobacteria
    ## Taxa_00158                         c__Actinobacteria
    ## Taxa_00159                         c__Actinobacteria
    ## Taxa_00160                         c__Actinobacteria
    ## Taxa_00161                         c__Actinobacteria
    ## Taxa_00162                         c__Actinobacteria
    ## Taxa_00163                         c__Actinobacteria
    ## Taxa_00164                         c__Actinobacteria
    ## Taxa_00165                         c__Actinobacteria
    ## Taxa_00166                         c__Actinobacteria
    ## Taxa_00167                         c__Actinobacteria
    ## Taxa_00168                         c__Actinobacteria
    ## Taxa_00169                         c__Actinobacteria
    ## Taxa_00170                         c__Actinobacteria
    ## Taxa_00171                         c__Actinobacteria
    ## Taxa_00172                         c__Actinobacteria
    ## Taxa_00173                         c__Actinobacteria
    ## Taxa_00174                         c__Actinobacteria
    ## Taxa_00175                         c__Actinobacteria
    ## Taxa_00176                         c__Actinobacteria
    ## Taxa_00177                         c__Actinobacteria
    ## Taxa_00178                         c__Actinobacteria
    ## Taxa_00179                         c__Actinobacteria
    ## Taxa_00180                         c__Actinobacteria
    ## Taxa_00181                         c__Actinobacteria
    ## Taxa_00182                         c__Actinobacteria
    ## Taxa_00183                         c__Actinobacteria
    ## Taxa_00184                         c__Actinobacteria
    ## Taxa_00185                         c__Actinobacteria
    ## Taxa_00186                         c__Actinobacteria
    ## Taxa_00187                         c__Actinobacteria
    ## Taxa_00188                         c__Actinobacteria
    ## Taxa_00189                         c__Actinobacteria
    ## Taxa_00190                         c__Actinobacteria
    ## Taxa_00191                         c__Actinobacteria
    ## Taxa_00192                         c__Actinobacteria
    ## Taxa_00193                         c__Actinobacteria
    ## Taxa_00194                         c__Actinobacteria
    ## Taxa_00195                         c__Actinobacteria
    ## Taxa_00196                         c__Actinobacteria
    ## Taxa_00197                         c__Actinobacteria
    ## Taxa_00198                         c__Actinobacteria
    ## Taxa_00199                         c__Actinobacteria
    ## Taxa_00200                         c__Actinobacteria
    ## Taxa_00201                         c__Actinobacteria
    ## Taxa_00202                         c__Actinobacteria
    ## Taxa_00203                         c__Actinobacteria
    ## Taxa_00204                         c__Actinobacteria
    ## Taxa_00205                         c__Actinobacteria
    ## Taxa_00206                         c__Actinobacteria
    ## Taxa_00207                         c__Actinobacteria
    ## Taxa_00208                         c__Actinobacteria
    ## Taxa_00209                         c__Actinobacteria
    ## Taxa_00210                         c__Actinobacteria
    ## Taxa_00211                         c__Actinobacteria
    ## Taxa_00212                         c__Actinobacteria
    ## Taxa_00213                         c__Actinobacteria
    ## Taxa_00214                         c__Actinobacteria
    ## Taxa_00215                         c__Actinobacteria
    ## Taxa_00216                         c__Actinobacteria
    ## Taxa_00217                        c__Thermoleophilia
    ## Taxa_00218                              c__Aquificae
    ## Taxa_00219                              c__Aquificae
    ## Taxa_00220                              c__Aquificae
    ## Taxa_00221                                      <NA>
    ## Taxa_00222                    c__Armatimonadetes_gp2
    ## Taxa_00223                    c__Armatimonadetes_gp4
    ## Taxa_00224                    c__Armatimonadetes_gp5
    ## Taxa_00225                          c__Armatimonadia
    ## Taxa_00226                       c__Chthonomonadetes
    ## Taxa_00227                         c__Fimbriimonadia
    ## Taxa_00228                                      <NA>
    ## Taxa_00229                            c__Bacteroidia
    ## Taxa_00230                            c__Bacteroidia
    ## Taxa_00231                            c__Bacteroidia
    ## Taxa_00232                            c__Bacteroidia
    ## Taxa_00233                            c__Bacteroidia
    ## Taxa_00234                             c__Cytophagia
    ## Taxa_00235                             c__Cytophagia
    ## Taxa_00236                             c__Cytophagia
    ## Taxa_00237                             c__Cytophagia
    ## Taxa_00238                             c__Cytophagia
    ## Taxa_00239                             c__Cytophagia
    ## Taxa_00240                             c__Cytophagia
    ## Taxa_00241                             c__Cytophagia
    ## Taxa_00242                             c__Cytophagia
    ## Taxa_00243                             c__Cytophagia
    ## Taxa_00244                             c__Cytophagia
    ## Taxa_00245                             c__Cytophagia
    ## Taxa_00246                             c__Cytophagia
    ## Taxa_00247                             c__Cytophagia
    ## Taxa_00248                             c__Cytophagia
    ## Taxa_00249                             c__Cytophagia
    ## Taxa_00250                             c__Cytophagia
    ## Taxa_00251                             c__Cytophagia
    ## Taxa_00252                             c__Cytophagia
    ## Taxa_00253                             c__Cytophagia
    ## Taxa_00254                             c__Cytophagia
    ## Taxa_00255                             c__Cytophagia
    ## Taxa_00256                             c__Cytophagia
    ## Taxa_00257                         c__Flavobacteriia
    ## Taxa_00258                         c__Flavobacteriia
    ## Taxa_00259                         c__Flavobacteriia
    ## Taxa_00260                         c__Flavobacteriia
    ## Taxa_00261                         c__Flavobacteriia
    ## Taxa_00262                         c__Flavobacteriia
    ## Taxa_00263                         c__Flavobacteriia
    ## Taxa_00264                         c__Flavobacteriia
    ## Taxa_00265                         c__Flavobacteriia
    ## Taxa_00266                         c__Flavobacteriia
    ## Taxa_00267                         c__Flavobacteriia
    ## Taxa_00268                         c__Flavobacteriia
    ## Taxa_00269                         c__Flavobacteriia
    ## Taxa_00270                       c__Sphingobacteriia
    ## Taxa_00271                       c__Sphingobacteriia
    ## Taxa_00272                       c__Sphingobacteriia
    ## Taxa_00273                       c__Sphingobacteriia
    ## Taxa_00274                       c__Sphingobacteriia
    ## Taxa_00275                       c__Sphingobacteriia
    ## Taxa_00276                       c__Sphingobacteriia
    ## Taxa_00277                       c__Sphingobacteriia
    ## Taxa_00278                       c__Sphingobacteriia
    ## Taxa_00279                       c__Sphingobacteriia
    ## Taxa_00280                       c__Sphingobacteriia
    ## Taxa_00281                       c__Sphingobacteriia
    ## Taxa_00282                       c__Sphingobacteriia
    ## Taxa_00283                       c__Sphingobacteriia
    ## Taxa_00284                       c__Sphingobacteriia
    ## Taxa_00285                       c__Sphingobacteriia
    ## Taxa_00286                       c__Sphingobacteriia
    ## Taxa_00287                       c__Sphingobacteriia
    ## Taxa_00288                       c__Sphingobacteriia
    ## Taxa_00289                       c__Sphingobacteriia
    ## Taxa_00290                       c__Sphingobacteriia
    ## Taxa_00291                       c__Sphingobacteriia
    ## Taxa_00292                       c__Sphingobacteriia
    ## Taxa_00293                       c__Sphingobacteriia
    ## Taxa_00294                       c__Sphingobacteriia
    ## Taxa_00295                       c__Sphingobacteriia
    ## Taxa_00296                       c__Sphingobacteriia
    ## Taxa_00297                       c__Sphingobacteriia
    ## Taxa_00298                       c__Sphingobacteriia
    ## Taxa_00299                       c__Sphingobacteriia
    ## Taxa_00300                       c__Sphingobacteriia
    ## Taxa_00301                       c__Sphingobacteriia
    ## Taxa_00302                       c__Sphingobacteriia
    ## Taxa_00303                       c__Sphingobacteriia
    ## Taxa_00304             c__BRC1_genera_incertae_sedis
    ## Taxa_00305            c__WPS-1_genera_incertae_sedis
    ## Taxa_00306            c__WPS-2_genera_incertae_sedis
    ## Taxa_00307              c__ZB3_genera_incertae_sedis
    ## Taxa_00308 c__Saccharibacteria_genera_incertae_sedis
    ## Taxa_00309                             c__Chlamydiia
    ## Taxa_00310                             c__Chlamydiia
    ## Taxa_00311                             c__Chlamydiia
    ## Taxa_00312                             c__Chlamydiia
    ## Taxa_00313                             c__Chlamydiia
    ## Taxa_00314                                      <NA>
    ## Taxa_00315                           c__Anaerolineae
    ## Taxa_00316                           c__Anaerolineae
    ## Taxa_00317                           c__Anaerolineae
    ## Taxa_00318                           c__Anaerolineae
    ## Taxa_00319                           c__Anaerolineae
    ## Taxa_00320                           c__Anaerolineae
    ## Taxa_00321                         c__Ardenticatenia
    ## Taxa_00322                            c__Caldilineae
    ## Taxa_00323                            c__Caldilineae
    ## Taxa_00324                            c__Caldilineae
    ## Taxa_00325                           c__Chloroflexia
    ## Taxa_00326                           c__Chloroflexia
    ## Taxa_00327                           c__Chloroflexia
    ## Taxa_00328                           c__Chloroflexia
    ## Taxa_00329                           c__Chloroflexia
    ## Taxa_00330                           c__Chloroflexia
    ## Taxa_00331                        c__Dehalococcoidia
    ## Taxa_00332                        c__Ktedonobacteria
    ## Taxa_00333                        c__Ktedonobacteria
    ## Taxa_00334                        c__Ktedonobacteria
    ## Taxa_00335                        c__Ktedonobacteria
    ## Taxa_00336                           c__Thermoflexia
    ## Taxa_00337                         c__Thermomicrobia
    ## Taxa_00338                         c__Thermomicrobia
    ## Taxa_00339                         c__Thermomicrobia
    ## Taxa_00340                         c__Thermomicrobia
    ## Taxa_00341                                      <NA>
    ## Taxa_00342                            c__Chloroplast
    ## Taxa_00343                            c__Chloroplast
    ## Taxa_00344                            c__Chloroplast
    ## Taxa_00345                            c__Chloroplast
    ## Taxa_00346                            c__Chloroplast
    ## Taxa_00347                            c__Chloroplast
    ## Taxa_00348                          c__Cyanobacteria
    ## Taxa_00349                             c__Deinococci
    ## Taxa_00350                                      <NA>
    ## Taxa_00351                          c__Elusimicrobia
    ## Taxa_00352                           c__Endomicrobia
    ## Taxa_00353                                      <NA>
    ## Taxa_00354                                c__Bacilli
    ## Taxa_00355                                c__Bacilli
    ## Taxa_00356                                c__Bacilli
    ## Taxa_00357                                c__Bacilli
    ## Taxa_00358                                c__Bacilli
    ## Taxa_00359                                c__Bacilli
    ## Taxa_00360                                c__Bacilli
    ## Taxa_00361                                c__Bacilli
    ## Taxa_00362                                c__Bacilli
    ## Taxa_00363                                c__Bacilli
    ## Taxa_00364                                c__Bacilli
    ## Taxa_00365                                c__Bacilli
    ## Taxa_00366                                c__Bacilli
    ## Taxa_00367                                c__Bacilli
    ## Taxa_00368                                c__Bacilli
    ## Taxa_00369                                c__Bacilli
    ## Taxa_00370                                c__Bacilli
    ## Taxa_00371                                c__Bacilli
    ## Taxa_00372                                c__Bacilli
    ## Taxa_00373                                c__Bacilli
    ## Taxa_00374                                c__Bacilli
    ## Taxa_00375                                c__Bacilli
    ## Taxa_00376                                c__Bacilli
    ## Taxa_00377                                c__Bacilli
    ## Taxa_00378                                c__Bacilli
    ## Taxa_00379                                c__Bacilli
    ## Taxa_00380                                c__Bacilli
    ## Taxa_00381                                c__Bacilli
    ## Taxa_00382                                c__Bacilli
    ## Taxa_00383                                c__Bacilli
    ## Taxa_00384                                c__Bacilli
    ## Taxa_00385                                c__Bacilli
    ## Taxa_00386                                c__Bacilli
    ## Taxa_00387                                c__Bacilli
    ## Taxa_00388                                c__Bacilli
    ## Taxa_00389                                c__Bacilli
    ## Taxa_00390                                c__Bacilli
    ## Taxa_00391                                c__Bacilli
    ## Taxa_00392                                c__Bacilli
    ## Taxa_00393                                c__Bacilli
    ## Taxa_00394                                c__Bacilli
    ## Taxa_00395                                c__Bacilli
    ## Taxa_00396                                c__Bacilli
    ## Taxa_00397                                c__Bacilli
    ## Taxa_00398                                c__Bacilli
    ## Taxa_00399                                c__Bacilli
    ## Taxa_00400                                c__Bacilli
    ## Taxa_00401                                c__Bacilli
    ## Taxa_00402                                c__Bacilli
    ## Taxa_00403                                c__Bacilli
    ## Taxa_00404                                c__Bacilli
    ## Taxa_00405                                c__Bacilli
    ## Taxa_00406                                c__Bacilli
    ## Taxa_00407                                c__Bacilli
    ## Taxa_00408                                c__Bacilli
    ## Taxa_00409                                c__Bacilli
    ## Taxa_00410                                c__Bacilli
    ## Taxa_00411                                c__Bacilli
    ## Taxa_00412                                c__Bacilli
    ## Taxa_00413                                c__Bacilli
    ## Taxa_00414                                c__Bacilli
    ## Taxa_00415                                c__Bacilli
    ## Taxa_00416                                c__Bacilli
    ## Taxa_00417                                c__Bacilli
    ## Taxa_00418                                c__Bacilli
    ## Taxa_00419                                c__Bacilli
    ## Taxa_00420                                c__Bacilli
    ## Taxa_00421                                c__Bacilli
    ## Taxa_00422                                c__Bacilli
    ## Taxa_00423                                c__Bacilli
    ## Taxa_00424                             c__Clostridia
    ## Taxa_00425                             c__Clostridia
    ## Taxa_00426                             c__Clostridia
    ## Taxa_00427                             c__Clostridia
    ## Taxa_00428                             c__Clostridia
    ## Taxa_00429                             c__Clostridia
    ## Taxa_00430                             c__Clostridia
    ## Taxa_00431                             c__Clostridia
    ## Taxa_00432                             c__Clostridia
    ## Taxa_00433                             c__Clostridia
    ## Taxa_00434                             c__Clostridia
    ## Taxa_00435                             c__Clostridia
    ## Taxa_00436                             c__Clostridia
    ## Taxa_00437                             c__Clostridia
    ## Taxa_00438                             c__Clostridia
    ## Taxa_00439                             c__Clostridia
    ## Taxa_00440                             c__Clostridia
    ## Taxa_00441                             c__Clostridia
    ## Taxa_00442                             c__Clostridia
    ## Taxa_00443                             c__Clostridia
    ## Taxa_00444                             c__Clostridia
    ## Taxa_00445                             c__Clostridia
    ## Taxa_00446                             c__Clostridia
    ## Taxa_00447                             c__Clostridia
    ## Taxa_00448                             c__Clostridia
    ## Taxa_00449                             c__Clostridia
    ## Taxa_00450                             c__Clostridia
    ## Taxa_00451                             c__Clostridia
    ## Taxa_00452                             c__Clostridia
    ## Taxa_00453                             c__Clostridia
    ## Taxa_00454                             c__Clostridia
    ## Taxa_00455                             c__Clostridia
    ## Taxa_00456                             c__Clostridia
    ## Taxa_00457                             c__Clostridia
    ## Taxa_00458                             c__Clostridia
    ## Taxa_00459                             c__Clostridia
    ## Taxa_00460                             c__Clostridia
    ## Taxa_00461                       c__Erysipelotrichia
    ## Taxa_00462                          c__Negativicutes
    ## Taxa_00463                          c__Negativicutes
    ## Taxa_00464                          c__Negativicutes
    ## Taxa_00465                          c__Negativicutes
    ## Taxa_00466                          c__Negativicutes
    ## Taxa_00467                       c__Gemmatimonadetes
    ## Taxa_00468               c__Candidatus Hydrogenedens
    ## Taxa_00469  c__Latescibacteria_genera_incertae_sedis
    ## Taxa_00470                                      <NA>
    ## Taxa_00471                          c__Lentisphaeria
    ## Taxa_00472   c__Microgenomates_genera_incertae_sedis
    ## Taxa_00473                             c__Nitrospira
    ## Taxa_00474     c__Omnitrophica_genera_incertae_sedis
    ## Taxa_00475    c__Parcubacteria_genera_incertae_sedis
    ## Taxa_00476                                      <NA>
    ## Taxa_00477                          c__Phycisphaerae
    ## Taxa_00478                          c__Phycisphaerae
    ## Taxa_00479                          c__Phycisphaerae
    ## Taxa_00480                         c__Planctomycetia
    ## Taxa_00481                         c__Planctomycetia
    ## Taxa_00482                         c__Planctomycetia
    ## Taxa_00483                         c__Planctomycetia
    ## Taxa_00484                         c__Planctomycetia
    ## Taxa_00485                         c__Planctomycetia
    ## Taxa_00486                         c__Planctomycetia
    ## Taxa_00487                         c__Planctomycetia
    ## Taxa_00488                         c__Planctomycetia
    ## Taxa_00489                         c__Planctomycetia
    ## Taxa_00490                         c__Planctomycetia
    ## Taxa_00491                         c__Planctomycetia
    ## Taxa_00492                         c__Planctomycetia
    ## Taxa_00493                         c__Planctomycetia
    ## Taxa_00494                         c__Planctomycetia
    ## Taxa_00495                         c__Planctomycetia
    ## Taxa_00496                         c__Planctomycetia
    ## Taxa_00497                         c__Planctomycetia
    ## Taxa_00498                                      <NA>
    ## Taxa_00499                    c__Alphaproteobacteria
    ## Taxa_00500                    c__Alphaproteobacteria
    ## Taxa_00501                    c__Alphaproteobacteria
    ## Taxa_00502                    c__Alphaproteobacteria
    ## Taxa_00503                    c__Alphaproteobacteria
    ## Taxa_00504                    c__Alphaproteobacteria
    ## Taxa_00505                    c__Alphaproteobacteria
    ## Taxa_00506                    c__Alphaproteobacteria
    ## Taxa_00507                    c__Alphaproteobacteria
    ## Taxa_00508                    c__Alphaproteobacteria
    ## Taxa_00509                    c__Alphaproteobacteria
    ## Taxa_00510                    c__Alphaproteobacteria
    ## Taxa_00511                    c__Alphaproteobacteria
    ## Taxa_00512                    c__Alphaproteobacteria
    ## Taxa_00513                    c__Alphaproteobacteria
    ## Taxa_00514                    c__Alphaproteobacteria
    ## Taxa_00515                    c__Alphaproteobacteria
    ## Taxa_00516                    c__Alphaproteobacteria
    ## Taxa_00517                    c__Alphaproteobacteria
    ## Taxa_00518                    c__Alphaproteobacteria
    ## Taxa_00519                    c__Alphaproteobacteria
    ## Taxa_00520                    c__Alphaproteobacteria
    ## Taxa_00521                    c__Alphaproteobacteria
    ## Taxa_00522                    c__Alphaproteobacteria
    ## Taxa_00523                    c__Alphaproteobacteria
    ## Taxa_00524                    c__Alphaproteobacteria
    ## Taxa_00525                    c__Alphaproteobacteria
    ## Taxa_00526                    c__Alphaproteobacteria
    ## Taxa_00527                    c__Alphaproteobacteria
    ## Taxa_00528                    c__Alphaproteobacteria
    ## Taxa_00529                    c__Alphaproteobacteria
    ## Taxa_00530                    c__Alphaproteobacteria
    ## Taxa_00531                    c__Alphaproteobacteria
    ## Taxa_00532                    c__Alphaproteobacteria
    ## Taxa_00533                    c__Alphaproteobacteria
    ## Taxa_00534                    c__Alphaproteobacteria
    ## Taxa_00535                    c__Alphaproteobacteria
    ## Taxa_00536                    c__Alphaproteobacteria
    ## Taxa_00537                    c__Alphaproteobacteria
    ## Taxa_00538                    c__Alphaproteobacteria
    ## Taxa_00539                    c__Alphaproteobacteria
    ## Taxa_00540                    c__Alphaproteobacteria
    ## Taxa_00541                    c__Alphaproteobacteria
    ## Taxa_00542                    c__Alphaproteobacteria
    ## Taxa_00543                    c__Alphaproteobacteria
    ## Taxa_00544                    c__Alphaproteobacteria
    ## Taxa_00545                    c__Alphaproteobacteria
    ## Taxa_00546                    c__Alphaproteobacteria
    ## Taxa_00547                    c__Alphaproteobacteria
    ## Taxa_00548                    c__Alphaproteobacteria
    ## Taxa_00549                    c__Alphaproteobacteria
    ## Taxa_00550                    c__Alphaproteobacteria
    ## Taxa_00551                    c__Alphaproteobacteria
    ## Taxa_00552                    c__Alphaproteobacteria
    ## Taxa_00553                    c__Alphaproteobacteria
    ## Taxa_00554                    c__Alphaproteobacteria
    ## Taxa_00555                    c__Alphaproteobacteria
    ## Taxa_00556                    c__Alphaproteobacteria
    ## Taxa_00557                    c__Alphaproteobacteria
    ## Taxa_00558                    c__Alphaproteobacteria
    ## Taxa_00559                    c__Alphaproteobacteria
    ## Taxa_00560                    c__Alphaproteobacteria
    ## Taxa_00561                    c__Alphaproteobacteria
    ## Taxa_00562                    c__Alphaproteobacteria
    ## Taxa_00563                    c__Alphaproteobacteria
    ## Taxa_00564                    c__Alphaproteobacteria
    ## Taxa_00565                    c__Alphaproteobacteria
    ## Taxa_00566                    c__Alphaproteobacteria
    ## Taxa_00567                    c__Alphaproteobacteria
    ## Taxa_00568                    c__Alphaproteobacteria
    ## Taxa_00569                    c__Alphaproteobacteria
    ## Taxa_00570                    c__Alphaproteobacteria
    ## Taxa_00571                    c__Alphaproteobacteria
    ## Taxa_00572                    c__Alphaproteobacteria
    ## Taxa_00573                    c__Alphaproteobacteria
    ## Taxa_00574                    c__Alphaproteobacteria
    ## Taxa_00575                    c__Alphaproteobacteria
    ## Taxa_00576                    c__Alphaproteobacteria
    ## Taxa_00577                    c__Alphaproteobacteria
    ## Taxa_00578                    c__Alphaproteobacteria
    ## Taxa_00579                    c__Alphaproteobacteria
    ## Taxa_00580                    c__Alphaproteobacteria
    ## Taxa_00581                    c__Alphaproteobacteria
    ## Taxa_00582                    c__Alphaproteobacteria
    ## Taxa_00583                    c__Alphaproteobacteria
    ## Taxa_00584                    c__Alphaproteobacteria
    ## Taxa_00585                    c__Alphaproteobacteria
    ## Taxa_00586                    c__Alphaproteobacteria
    ## Taxa_00587                    c__Alphaproteobacteria
    ## Taxa_00588                    c__Alphaproteobacteria
    ## Taxa_00589                    c__Alphaproteobacteria
    ## Taxa_00590                    c__Alphaproteobacteria
    ## Taxa_00591                    c__Alphaproteobacteria
    ## Taxa_00592                    c__Alphaproteobacteria
    ## Taxa_00593                    c__Alphaproteobacteria
    ## Taxa_00594                    c__Alphaproteobacteria
    ## Taxa_00595                    c__Alphaproteobacteria
    ## Taxa_00596                    c__Alphaproteobacteria
    ## Taxa_00597                    c__Alphaproteobacteria
    ## Taxa_00598                    c__Alphaproteobacteria
    ## Taxa_00599                    c__Alphaproteobacteria
    ## Taxa_00600                    c__Alphaproteobacteria
    ## Taxa_00601                    c__Alphaproteobacteria
    ## Taxa_00602                    c__Alphaproteobacteria
    ## Taxa_00603                    c__Alphaproteobacteria
    ## Taxa_00604                    c__Alphaproteobacteria
    ## Taxa_00605                    c__Alphaproteobacteria
    ## Taxa_00606                    c__Alphaproteobacteria
    ## Taxa_00607                    c__Alphaproteobacteria
    ## Taxa_00608                    c__Alphaproteobacteria
    ## Taxa_00609                    c__Alphaproteobacteria
    ## Taxa_00610                    c__Alphaproteobacteria
    ## Taxa_00611                    c__Alphaproteobacteria
    ## Taxa_00612                    c__Alphaproteobacteria
    ## Taxa_00613                    c__Alphaproteobacteria
    ## Taxa_00614                    c__Alphaproteobacteria
    ## Taxa_00615                    c__Alphaproteobacteria
    ## Taxa_00616                    c__Alphaproteobacteria
    ## Taxa_00617                    c__Alphaproteobacteria
    ## Taxa_00618                    c__Alphaproteobacteria
    ## Taxa_00619                    c__Alphaproteobacteria
    ## Taxa_00620                    c__Alphaproteobacteria
    ## Taxa_00621                    c__Alphaproteobacteria
    ## Taxa_00622                    c__Alphaproteobacteria
    ## Taxa_00623                    c__Alphaproteobacteria
    ## Taxa_00624                    c__Alphaproteobacteria
    ## Taxa_00625                    c__Alphaproteobacteria
    ## Taxa_00626                    c__Alphaproteobacteria
    ## Taxa_00627                    c__Alphaproteobacteria
    ## Taxa_00628                    c__Alphaproteobacteria
    ## Taxa_00629                    c__Alphaproteobacteria
    ## Taxa_00630                    c__Alphaproteobacteria
    ## Taxa_00631                    c__Alphaproteobacteria
    ## Taxa_00632                    c__Alphaproteobacteria
    ## Taxa_00633                    c__Alphaproteobacteria
    ## Taxa_00634                    c__Alphaproteobacteria
    ## Taxa_00635                    c__Alphaproteobacteria
    ## Taxa_00636                    c__Alphaproteobacteria
    ## Taxa_00637                    c__Alphaproteobacteria
    ## Taxa_00638                    c__Alphaproteobacteria
    ## Taxa_00639                    c__Alphaproteobacteria
    ## Taxa_00640                    c__Alphaproteobacteria
    ## Taxa_00641                    c__Alphaproteobacteria
    ## Taxa_00642                    c__Alphaproteobacteria
    ## Taxa_00643                    c__Alphaproteobacteria
    ## Taxa_00644                    c__Alphaproteobacteria
    ## Taxa_00645                    c__Alphaproteobacteria
    ## Taxa_00646                    c__Alphaproteobacteria
    ## Taxa_00647                    c__Alphaproteobacteria
    ## Taxa_00648                    c__Alphaproteobacteria
    ## Taxa_00649                    c__Alphaproteobacteria
    ## Taxa_00650                    c__Alphaproteobacteria
    ## Taxa_00651                    c__Alphaproteobacteria
    ## Taxa_00652                    c__Alphaproteobacteria
    ## Taxa_00653                    c__Alphaproteobacteria
    ## Taxa_00654                    c__Alphaproteobacteria
    ## Taxa_00655                    c__Alphaproteobacteria
    ## Taxa_00656                    c__Alphaproteobacteria
    ## Taxa_00657                    c__Alphaproteobacteria
    ## Taxa_00658                    c__Alphaproteobacteria
    ## Taxa_00659                    c__Alphaproteobacteria
    ## Taxa_00660                    c__Alphaproteobacteria
    ## Taxa_00661                    c__Alphaproteobacteria
    ## Taxa_00662                    c__Alphaproteobacteria
    ## Taxa_00663                    c__Alphaproteobacteria
    ## Taxa_00664                    c__Alphaproteobacteria
    ## Taxa_00665                    c__Alphaproteobacteria
    ## Taxa_00666                    c__Alphaproteobacteria
    ## Taxa_00667                    c__Alphaproteobacteria
    ## Taxa_00668                    c__Alphaproteobacteria
    ## Taxa_00669                    c__Alphaproteobacteria
    ## Taxa_00670                    c__Alphaproteobacteria
    ## Taxa_00671                    c__Alphaproteobacteria
    ## Taxa_00672                    c__Alphaproteobacteria
    ## Taxa_00673                    c__Alphaproteobacteria
    ## Taxa_00674                    c__Alphaproteobacteria
    ## Taxa_00675                    c__Alphaproteobacteria
    ## Taxa_00676                    c__Alphaproteobacteria
    ## Taxa_00677                    c__Alphaproteobacteria
    ## Taxa_00678                    c__Alphaproteobacteria
    ## Taxa_00679                    c__Alphaproteobacteria
    ## Taxa_00680                    c__Alphaproteobacteria
    ## Taxa_00681                    c__Alphaproteobacteria
    ## Taxa_00682                     c__Betaproteobacteria
    ## Taxa_00683                     c__Betaproteobacteria
    ## Taxa_00684                     c__Betaproteobacteria
    ## Taxa_00685                     c__Betaproteobacteria
    ## Taxa_00686                     c__Betaproteobacteria
    ## Taxa_00687                     c__Betaproteobacteria
    ## Taxa_00688                     c__Betaproteobacteria
    ## Taxa_00689                     c__Betaproteobacteria
    ## Taxa_00690                     c__Betaproteobacteria
    ## Taxa_00691                     c__Betaproteobacteria
    ## Taxa_00692                     c__Betaproteobacteria
    ## Taxa_00693                     c__Betaproteobacteria
    ## Taxa_00694                     c__Betaproteobacteria
    ## Taxa_00695                     c__Betaproteobacteria
    ## Taxa_00696                     c__Betaproteobacteria
    ## Taxa_00697                     c__Betaproteobacteria
    ## Taxa_00698                     c__Betaproteobacteria
    ## Taxa_00699                     c__Betaproteobacteria
    ## Taxa_00700                     c__Betaproteobacteria
    ## Taxa_00701                     c__Betaproteobacteria
    ## Taxa_00702                     c__Betaproteobacteria
    ## Taxa_00703                     c__Betaproteobacteria
    ## Taxa_00704                     c__Betaproteobacteria
    ## Taxa_00705                     c__Betaproteobacteria
    ## Taxa_00706                     c__Betaproteobacteria
    ## Taxa_00707                     c__Betaproteobacteria
    ## Taxa_00708                     c__Betaproteobacteria
    ## Taxa_00709                     c__Betaproteobacteria
    ## Taxa_00710                     c__Betaproteobacteria
    ## Taxa_00711                     c__Betaproteobacteria
    ## Taxa_00712                     c__Betaproteobacteria
    ## Taxa_00713                     c__Betaproteobacteria
    ## Taxa_00714                     c__Betaproteobacteria
    ## Taxa_00715                     c__Betaproteobacteria
    ## Taxa_00716                     c__Betaproteobacteria
    ## Taxa_00717                     c__Betaproteobacteria
    ## Taxa_00718                     c__Betaproteobacteria
    ## Taxa_00719                     c__Betaproteobacteria
    ## Taxa_00720                     c__Betaproteobacteria
    ## Taxa_00721                     c__Betaproteobacteria
    ## Taxa_00722                     c__Betaproteobacteria
    ## Taxa_00723                     c__Betaproteobacteria
    ## Taxa_00724                     c__Betaproteobacteria
    ## Taxa_00725                     c__Betaproteobacteria
    ## Taxa_00726                     c__Betaproteobacteria
    ## Taxa_00727                     c__Betaproteobacteria
    ## Taxa_00728                     c__Betaproteobacteria
    ## Taxa_00729                     c__Betaproteobacteria
    ## Taxa_00730                     c__Betaproteobacteria
    ## Taxa_00731                     c__Betaproteobacteria
    ## Taxa_00732                     c__Betaproteobacteria
    ## Taxa_00733                     c__Betaproteobacteria
    ## Taxa_00734                     c__Betaproteobacteria
    ## Taxa_00735                     c__Betaproteobacteria
    ## Taxa_00736                     c__Betaproteobacteria
    ## Taxa_00737                     c__Betaproteobacteria
    ## Taxa_00738                     c__Betaproteobacteria
    ## Taxa_00739                     c__Betaproteobacteria
    ## Taxa_00740                     c__Betaproteobacteria
    ## Taxa_00741                     c__Betaproteobacteria
    ## Taxa_00742                     c__Betaproteobacteria
    ## Taxa_00743                     c__Betaproteobacteria
    ## Taxa_00744                     c__Betaproteobacteria
    ## Taxa_00745                     c__Betaproteobacteria
    ## Taxa_00746                     c__Betaproteobacteria
    ## Taxa_00747                     c__Betaproteobacteria
    ## Taxa_00748                     c__Betaproteobacteria
    ## Taxa_00749                     c__Betaproteobacteria
    ## Taxa_00750                     c__Betaproteobacteria
    ## Taxa_00751                     c__Betaproteobacteria
    ## Taxa_00752                     c__Betaproteobacteria
    ## Taxa_00753                     c__Betaproteobacteria
    ## Taxa_00754                     c__Betaproteobacteria
    ## Taxa_00755                     c__Betaproteobacteria
    ## Taxa_00756                     c__Betaproteobacteria
    ## Taxa_00757                     c__Betaproteobacteria
    ## Taxa_00758                     c__Betaproteobacteria
    ## Taxa_00759                     c__Betaproteobacteria
    ## Taxa_00760                     c__Betaproteobacteria
    ## Taxa_00761                     c__Betaproteobacteria
    ## Taxa_00762                     c__Betaproteobacteria
    ## Taxa_00763                     c__Betaproteobacteria
    ## Taxa_00764                     c__Betaproteobacteria
    ## Taxa_00765                     c__Betaproteobacteria
    ## Taxa_00766                     c__Betaproteobacteria
    ## Taxa_00767                     c__Betaproteobacteria
    ## Taxa_00768                     c__Betaproteobacteria
    ## Taxa_00769                     c__Betaproteobacteria
    ## Taxa_00770                     c__Betaproteobacteria
    ## Taxa_00771                     c__Betaproteobacteria
    ## Taxa_00772                     c__Betaproteobacteria
    ## Taxa_00773                     c__Betaproteobacteria
    ## Taxa_00774                     c__Betaproteobacteria
    ## Taxa_00775                     c__Betaproteobacteria
    ## Taxa_00776                     c__Betaproteobacteria
    ## Taxa_00777                     c__Betaproteobacteria
    ## Taxa_00778                     c__Betaproteobacteria
    ## Taxa_00779                    c__Deltaproteobacteria
    ## Taxa_00780                    c__Deltaproteobacteria
    ## Taxa_00781                    c__Deltaproteobacteria
    ## Taxa_00782                    c__Deltaproteobacteria
    ## Taxa_00783                    c__Deltaproteobacteria
    ## Taxa_00784                    c__Deltaproteobacteria
    ## Taxa_00785                    c__Deltaproteobacteria
    ## Taxa_00786                    c__Deltaproteobacteria
    ## Taxa_00787                    c__Deltaproteobacteria
    ## Taxa_00788                    c__Deltaproteobacteria
    ## Taxa_00789                    c__Deltaproteobacteria
    ## Taxa_00790                    c__Deltaproteobacteria
    ## Taxa_00791                    c__Deltaproteobacteria
    ## Taxa_00792                    c__Deltaproteobacteria
    ## Taxa_00793                    c__Deltaproteobacteria
    ## Taxa_00794                    c__Deltaproteobacteria
    ## Taxa_00795                    c__Deltaproteobacteria
    ## Taxa_00796                    c__Deltaproteobacteria
    ## Taxa_00797                    c__Deltaproteobacteria
    ## Taxa_00798                    c__Deltaproteobacteria
    ## Taxa_00799                    c__Deltaproteobacteria
    ## Taxa_00800                    c__Deltaproteobacteria
    ## Taxa_00801                    c__Deltaproteobacteria
    ## Taxa_00802                    c__Deltaproteobacteria
    ## Taxa_00803                    c__Deltaproteobacteria
    ## Taxa_00804                    c__Deltaproteobacteria
    ## Taxa_00805                    c__Deltaproteobacteria
    ## Taxa_00806                    c__Deltaproteobacteria
    ## Taxa_00807                    c__Deltaproteobacteria
    ## Taxa_00808                    c__Deltaproteobacteria
    ## Taxa_00809                    c__Deltaproteobacteria
    ## Taxa_00810                    c__Deltaproteobacteria
    ## Taxa_00811                    c__Deltaproteobacteria
    ## Taxa_00812                    c__Deltaproteobacteria
    ## Taxa_00813                    c__Deltaproteobacteria
    ## Taxa_00814                    c__Deltaproteobacteria
    ## Taxa_00815                    c__Deltaproteobacteria
    ## Taxa_00816                    c__Deltaproteobacteria
    ## Taxa_00817                    c__Deltaproteobacteria
    ## Taxa_00818                    c__Deltaproteobacteria
    ## Taxa_00819                    c__Deltaproteobacteria
    ## Taxa_00820                    c__Deltaproteobacteria
    ## Taxa_00821                    c__Deltaproteobacteria
    ## Taxa_00822                    c__Deltaproteobacteria
    ## Taxa_00823                    c__Deltaproteobacteria
    ## Taxa_00824                    c__Deltaproteobacteria
    ## Taxa_00825                    c__Deltaproteobacteria
    ## Taxa_00826                    c__Deltaproteobacteria
    ## Taxa_00827                    c__Deltaproteobacteria
    ## Taxa_00828                    c__Deltaproteobacteria
    ## Taxa_00829                    c__Deltaproteobacteria
    ## Taxa_00830                  c__Epsilonproteobacteria
    ## Taxa_00831                  c__Epsilonproteobacteria
    ## Taxa_00832                    c__Gammaproteobacteria
    ## Taxa_00833                    c__Gammaproteobacteria
    ## Taxa_00834                    c__Gammaproteobacteria
    ## Taxa_00835                    c__Gammaproteobacteria
    ## Taxa_00836                    c__Gammaproteobacteria
    ## Taxa_00837                    c__Gammaproteobacteria
    ## Taxa_00838                    c__Gammaproteobacteria
    ## Taxa_00839                    c__Gammaproteobacteria
    ## Taxa_00840                    c__Gammaproteobacteria
    ## Taxa_00841                    c__Gammaproteobacteria
    ## Taxa_00842                    c__Gammaproteobacteria
    ## Taxa_00843                    c__Gammaproteobacteria
    ## Taxa_00844                    c__Gammaproteobacteria
    ## Taxa_00845                    c__Gammaproteobacteria
    ## Taxa_00846                    c__Gammaproteobacteria
    ## Taxa_00847                    c__Gammaproteobacteria
    ## Taxa_00848                    c__Gammaproteobacteria
    ## Taxa_00849                    c__Gammaproteobacteria
    ## Taxa_00850                    c__Gammaproteobacteria
    ## Taxa_00851                    c__Gammaproteobacteria
    ## Taxa_00852                    c__Gammaproteobacteria
    ## Taxa_00853                    c__Gammaproteobacteria
    ## Taxa_00854                    c__Gammaproteobacteria
    ## Taxa_00855                    c__Gammaproteobacteria
    ## Taxa_00856                    c__Gammaproteobacteria
    ## Taxa_00857                    c__Gammaproteobacteria
    ## Taxa_00858                    c__Gammaproteobacteria
    ## Taxa_00859                    c__Gammaproteobacteria
    ## Taxa_00860                    c__Gammaproteobacteria
    ## Taxa_00861                    c__Gammaproteobacteria
    ## Taxa_00862                    c__Gammaproteobacteria
    ## Taxa_00863                    c__Gammaproteobacteria
    ## Taxa_00864                    c__Gammaproteobacteria
    ## Taxa_00865                    c__Gammaproteobacteria
    ## Taxa_00866                    c__Gammaproteobacteria
    ## Taxa_00867                    c__Gammaproteobacteria
    ## Taxa_00868                    c__Gammaproteobacteria
    ## Taxa_00869                    c__Gammaproteobacteria
    ## Taxa_00870                    c__Gammaproteobacteria
    ## Taxa_00871                    c__Gammaproteobacteria
    ## Taxa_00872                    c__Gammaproteobacteria
    ## Taxa_00873                    c__Gammaproteobacteria
    ## Taxa_00874                    c__Gammaproteobacteria
    ## Taxa_00875                    c__Gammaproteobacteria
    ## Taxa_00876                    c__Gammaproteobacteria
    ## Taxa_00877                    c__Gammaproteobacteria
    ## Taxa_00878                    c__Gammaproteobacteria
    ## Taxa_00879                    c__Gammaproteobacteria
    ## Taxa_00880                    c__Gammaproteobacteria
    ## Taxa_00881                    c__Gammaproteobacteria
    ## Taxa_00882                    c__Gammaproteobacteria
    ## Taxa_00883                    c__Gammaproteobacteria
    ## Taxa_00884                    c__Gammaproteobacteria
    ## Taxa_00885                    c__Gammaproteobacteria
    ## Taxa_00886                    c__Gammaproteobacteria
    ## Taxa_00887                    c__Gammaproteobacteria
    ## Taxa_00888                    c__Gammaproteobacteria
    ## Taxa_00889                    c__Gammaproteobacteria
    ## Taxa_00890                    c__Gammaproteobacteria
    ## Taxa_00891                    c__Gammaproteobacteria
    ## Taxa_00892                    c__Gammaproteobacteria
    ## Taxa_00893                    c__Gammaproteobacteria
    ## Taxa_00894                    c__Gammaproteobacteria
    ## Taxa_00895                    c__Gammaproteobacteria
    ## Taxa_00896                    c__Gammaproteobacteria
    ## Taxa_00897                    c__Gammaproteobacteria
    ## Taxa_00898                    c__Gammaproteobacteria
    ## Taxa_00899                    c__Gammaproteobacteria
    ## Taxa_00900                    c__Gammaproteobacteria
    ## Taxa_00901                    c__Gammaproteobacteria
    ## Taxa_00902                    c__Gammaproteobacteria
    ## Taxa_00903                    c__Gammaproteobacteria
    ## Taxa_00904                    c__Gammaproteobacteria
    ## Taxa_00905                    c__Gammaproteobacteria
    ## Taxa_00906                    c__Gammaproteobacteria
    ## Taxa_00907                    c__Gammaproteobacteria
    ## Taxa_00908                    c__Gammaproteobacteria
    ## Taxa_00909                    c__Gammaproteobacteria
    ## Taxa_00910                    c__Gammaproteobacteria
    ## Taxa_00911                    c__Gammaproteobacteria
    ## Taxa_00912                    c__Gammaproteobacteria
    ## Taxa_00913                    c__Gammaproteobacteria
    ## Taxa_00914                            c__Oligoflexia
    ## Taxa_00915                           c__Spirochaetia
    ## Taxa_00916                           c__Spirochaetia
    ## Taxa_00917                           c__Spirochaetia
    ## Taxa_00918                           c__Spirochaetia
    ## Taxa_00919                           c__Spirochaetia
    ## Taxa_00920                           c__Spirochaetia
    ## Taxa_00921                           c__Spirochaetia
    ## Taxa_00922                           c__Spirochaetia
    ## Taxa_00923              c__SR1_genera_incertae_sedis
    ## Taxa_00924                             c__Mollicutes
    ## Taxa_00925                             c__Mollicutes
    ## Taxa_00926                  c__Thermodesulfobacteria
    ## Taxa_00927                                      <NA>
    ## Taxa_00928                               c__Opitutae
    ## Taxa_00929                               c__Opitutae
    ## Taxa_00930                               c__Opitutae
    ## Taxa_00931                               c__Opitutae
    ## Taxa_00932                               c__Opitutae
    ## Taxa_00933                               c__Opitutae
    ## Taxa_00934                               c__Opitutae
    ## Taxa_00935                         c__Spartobacteria
    ## Taxa_00936                         c__Spartobacteria
    ## Taxa_00937                         c__Spartobacteria
    ## Taxa_00938                           c__Subdivision3
    ## Taxa_00939                           c__Subdivision3
    ## Taxa_00940                           c__Subdivision3
    ## Taxa_00941                       c__Verrucomicrobiae
    ## Taxa_00942                       c__Verrucomicrobiae
    ## Taxa_00943                       c__Verrucomicrobiae
    ## Taxa_00944                       c__Verrucomicrobiae
    ## Taxa_00945                       c__Verrucomicrobiae
    ## Taxa_00946                       c__Verrucomicrobiae
    ## Taxa_00947                       c__Verrucomicrobiae
    ## Taxa_00948                       c__Verrucomicrobiae
    ## Taxa_00949                                      <NA>
    ##                                                Order
    ## Taxa_00000                                      <NA>
    ## Taxa_00001                      o__Desulfurococcales
    ## Taxa_00002                                      <NA>
    ## Taxa_00003                                      <NA>
    ## Taxa_00004                                      <NA>
    ## Taxa_00005                                   o__Gp10
    ## Taxa_00006                                   o__Gp11
    ## Taxa_00007                                   o__Gp12
    ## Taxa_00008                                   o__Gp13
    ## Taxa_00009                                   o__Gp15
    ## Taxa_00010                                   o__Gp16
    ## Taxa_00011                                   o__Gp17
    ## Taxa_00012                                   o__Gp18
    ## Taxa_00013                                   o__Gp19
    ## Taxa_00014                             o__Acidicapsa
    ## Taxa_00015                              o__Acidipila
    ## Taxa_00016                         o__Acidobacterium
    ## Taxa_00017                              o__Bryocella
    ## Taxa_00018                  o__Candidatus Koribacter
    ## Taxa_00019                           o__Edaphobacter
    ## Taxa_00020                                    o__Gp1
    ## Taxa_00021                           o__Granulicella
    ## Taxa_00022                          o__Telmatobacter
    ## Taxa_00023                            o__Terriglobus
    ## Taxa_00024                                   o__Gp20
    ## Taxa_00025                                   o__Gp22
    ## Taxa_00026                                   o__Gp25
    ## Taxa_00027                                    o__Gp2
    ## Taxa_00028                                      <NA>
    ## Taxa_00029                             o__Bryobacter
    ## Taxa_00030                  o__Candidatus Solibacter
    ## Taxa_00031                                    o__Gp3
    ## Taxa_00032                          o__Paludibaculum
    ## Taxa_00033                                      <NA>
    ## Taxa_00034                            o__Aridibacter
    ## Taxa_00035                          o__Blastocatella
    ## Taxa_00036                                    o__Gp4
    ## Taxa_00037                            o__Pyrinomonas
    ## Taxa_00038                                    o__Gp5
    ## Taxa_00039                                    o__Gp6
    ## Taxa_00040                                    o__Gp7
    ## Taxa_00041                           o__Holophagales
    ## Taxa_00042                           o__Holophagales
    ## Taxa_00043                           o__Holophagales
    ## Taxa_00044                                      <NA>
    ## Taxa_00045                                      <NA>
    ## Taxa_00046                       o__Acidimicrobiales
    ## Taxa_00047                       o__Acidimicrobiales
    ## Taxa_00048                       o__Acidimicrobiales
    ## Taxa_00049                       o__Acidimicrobiales
    ## Taxa_00050                       o__Acidimicrobiales
    ## Taxa_00051                       o__Acidimicrobiales
    ## Taxa_00052                       o__Acidimicrobiales
    ## Taxa_00053                       o__Acidimicrobiales
    ## Taxa_00054                       o__Acidimicrobiales
    ## Taxa_00055                        o__Actinomycetales
    ## Taxa_00056                        o__Actinomycetales
    ## Taxa_00057                        o__Actinomycetales
    ## Taxa_00058                        o__Actinomycetales
    ## Taxa_00059                        o__Actinomycetales
    ## Taxa_00060                        o__Actinomycetales
    ## Taxa_00061                        o__Actinomycetales
    ## Taxa_00062                        o__Actinomycetales
    ## Taxa_00063                        o__Actinomycetales
    ## Taxa_00064                        o__Actinomycetales
    ## Taxa_00065                        o__Actinomycetales
    ## Taxa_00066                        o__Actinomycetales
    ## Taxa_00067                        o__Actinomycetales
    ## Taxa_00068                        o__Actinomycetales
    ## Taxa_00069                        o__Actinomycetales
    ## Taxa_00070                        o__Actinomycetales
    ## Taxa_00071                        o__Actinomycetales
    ## Taxa_00072                        o__Actinomycetales
    ## Taxa_00073                        o__Actinomycetales
    ## Taxa_00074                        o__Actinomycetales
    ## Taxa_00075                        o__Actinomycetales
    ## Taxa_00076                        o__Actinomycetales
    ## Taxa_00077                        o__Actinomycetales
    ## Taxa_00078                        o__Actinomycetales
    ## Taxa_00079                        o__Actinomycetales
    ## Taxa_00080                        o__Actinomycetales
    ## Taxa_00081                        o__Actinomycetales
    ## Taxa_00082                        o__Actinomycetales
    ## Taxa_00083                        o__Actinomycetales
    ## Taxa_00084                        o__Actinomycetales
    ## Taxa_00085                        o__Actinomycetales
    ## Taxa_00086                        o__Actinomycetales
    ## Taxa_00087                        o__Actinomycetales
    ## Taxa_00088                        o__Actinomycetales
    ## Taxa_00089                        o__Actinomycetales
    ## Taxa_00090                        o__Actinomycetales
    ## Taxa_00091                        o__Actinomycetales
    ## Taxa_00092                        o__Actinomycetales
    ## Taxa_00093                        o__Actinomycetales
    ## Taxa_00094                        o__Actinomycetales
    ## Taxa_00095                        o__Actinomycetales
    ## Taxa_00096                        o__Actinomycetales
    ## Taxa_00097                        o__Actinomycetales
    ## Taxa_00098                        o__Actinomycetales
    ## Taxa_00099                        o__Actinomycetales
    ## Taxa_00100                        o__Actinomycetales
    ## Taxa_00101                        o__Actinomycetales
    ## Taxa_00102                        o__Actinomycetales
    ## Taxa_00103                        o__Actinomycetales
    ## Taxa_00104                        o__Actinomycetales
    ## Taxa_00105                        o__Actinomycetales
    ## Taxa_00106                        o__Actinomycetales
    ## Taxa_00107                        o__Actinomycetales
    ## Taxa_00108                        o__Actinomycetales
    ## Taxa_00109                        o__Actinomycetales
    ## Taxa_00110                        o__Actinomycetales
    ## Taxa_00111                        o__Actinomycetales
    ## Taxa_00112                        o__Actinomycetales
    ## Taxa_00113                        o__Actinomycetales
    ## Taxa_00114                        o__Actinomycetales
    ## Taxa_00115                        o__Actinomycetales
    ## Taxa_00116                        o__Actinomycetales
    ## Taxa_00117                        o__Actinomycetales
    ## Taxa_00118                        o__Actinomycetales
    ## Taxa_00119                        o__Actinomycetales
    ## Taxa_00120                        o__Actinomycetales
    ## Taxa_00121                        o__Actinomycetales
    ## Taxa_00122                        o__Actinomycetales
    ## Taxa_00123                        o__Actinomycetales
    ## Taxa_00124                        o__Actinomycetales
    ## Taxa_00125                        o__Actinomycetales
    ## Taxa_00126                        o__Actinomycetales
    ## Taxa_00127                        o__Actinomycetales
    ## Taxa_00128                        o__Actinomycetales
    ## Taxa_00129                        o__Actinomycetales
    ## Taxa_00130                        o__Actinomycetales
    ## Taxa_00131                        o__Actinomycetales
    ## Taxa_00132                        o__Actinomycetales
    ## Taxa_00133                        o__Actinomycetales
    ## Taxa_00134                        o__Actinomycetales
    ## Taxa_00135                        o__Actinomycetales
    ## Taxa_00136                        o__Actinomycetales
    ## Taxa_00137                        o__Actinomycetales
    ## Taxa_00138                        o__Actinomycetales
    ## Taxa_00139                        o__Actinomycetales
    ## Taxa_00140                        o__Actinomycetales
    ## Taxa_00141                        o__Actinomycetales
    ## Taxa_00142                        o__Actinomycetales
    ## Taxa_00143                        o__Actinomycetales
    ## Taxa_00144                        o__Actinomycetales
    ## Taxa_00145                        o__Actinomycetales
    ## Taxa_00146                        o__Actinomycetales
    ## Taxa_00147                        o__Actinomycetales
    ## Taxa_00148                        o__Actinomycetales
    ## Taxa_00149                        o__Actinomycetales
    ## Taxa_00150                        o__Actinomycetales
    ## Taxa_00151                        o__Actinomycetales
    ## Taxa_00152                        o__Actinomycetales
    ## Taxa_00153                        o__Actinomycetales
    ## Taxa_00154                        o__Actinomycetales
    ## Taxa_00155                        o__Actinomycetales
    ## Taxa_00156                        o__Actinomycetales
    ## Taxa_00157                        o__Actinomycetales
    ## Taxa_00158                        o__Actinomycetales
    ## Taxa_00159                        o__Actinomycetales
    ## Taxa_00160                        o__Actinomycetales
    ## Taxa_00161                        o__Actinomycetales
    ## Taxa_00162                        o__Actinomycetales
    ## Taxa_00163                        o__Actinomycetales
    ## Taxa_00164                        o__Actinomycetales
    ## Taxa_00165                        o__Actinomycetales
    ## Taxa_00166                        o__Actinomycetales
    ## Taxa_00167                        o__Actinomycetales
    ## Taxa_00168                        o__Actinomycetales
    ## Taxa_00169                        o__Actinomycetales
    ## Taxa_00170                        o__Actinomycetales
    ## Taxa_00171                        o__Actinomycetales
    ## Taxa_00172                        o__Actinomycetales
    ## Taxa_00173                        o__Actinomycetales
    ## Taxa_00174                        o__Actinomycetales
    ## Taxa_00175                        o__Actinomycetales
    ## Taxa_00176                        o__Actinomycetales
    ## Taxa_00177                        o__Actinomycetales
    ## Taxa_00178                        o__Actinomycetales
    ## Taxa_00179                        o__Actinomycetales
    ## Taxa_00180                        o__Actinomycetales
    ## Taxa_00181                        o__Actinomycetales
    ## Taxa_00182                        o__Actinomycetales
    ## Taxa_00183                        o__Actinomycetales
    ## Taxa_00184                        o__Actinomycetales
    ## Taxa_00185                        o__Actinomycetales
    ## Taxa_00186                        o__Actinomycetales
    ## Taxa_00187                        o__Actinomycetales
    ## Taxa_00188                        o__Actinomycetales
    ## Taxa_00189                        o__Actinomycetales
    ## Taxa_00190                        o__Actinomycetales
    ## Taxa_00191                        o__Actinomycetales
    ## Taxa_00192                        o__Actinomycetales
    ## Taxa_00193                        o__Actinomycetales
    ## Taxa_00194                        o__Actinomycetales
    ## Taxa_00195                        o__Actinomycetales
    ## Taxa_00196                        o__Actinomycetales
    ## Taxa_00197                        o__Actinomycetales
    ## Taxa_00198                        o__Actinomycetales
    ## Taxa_00199                        o__Actinomycetales
    ## Taxa_00200                        o__Actinomycetales
    ## Taxa_00201                        o__Actinomycetales
    ## Taxa_00202                        o__Actinomycetales
    ## Taxa_00203                        o__Actinomycetales
    ## Taxa_00204                        o__Actinomycetales
    ## Taxa_00205                        o__Actinomycetales
    ## Taxa_00206                        o__Actinomycetales
    ## Taxa_00207                        o__Actinomycetales
    ## Taxa_00208                        o__Actinomycetales
    ## Taxa_00209                        o__Actinomycetales
    ## Taxa_00210                        o__Actinomycetales
    ## Taxa_00211                        o__Actinomycetales
    ## Taxa_00212                             o__Gaiellales
    ## Taxa_00213                    o__Solirubrobacterales
    ## Taxa_00214                    o__Solirubrobacterales
    ## Taxa_00215                    o__Solirubrobacterales
    ## Taxa_00216                    o__Solirubrobacterales
    ## Taxa_00217                      o__Thermoleophilales
    ## Taxa_00218                                      <NA>
    ## Taxa_00219                            o__Aquificales
    ## Taxa_00220                            o__Aquificales
    ## Taxa_00221                                      <NA>
    ## Taxa_00222                    o__Armatimonadetes_gp2
    ## Taxa_00223                    o__Armatimonadetes_gp4
    ## Taxa_00224                    o__Armatimonadetes_gp5
    ## Taxa_00225                        o__Armatimonadales
    ## Taxa_00226                       o__Chthonomonadales
    ## Taxa_00227                       o__Fimbriimonadales
    ## Taxa_00228                                      <NA>
    ## Taxa_00229                          o__Bacteroidales
    ## Taxa_00230                          o__Bacteroidales
    ## Taxa_00231                          o__Bacteroidales
    ## Taxa_00232                          o__Bacteroidales
    ## Taxa_00233                          o__Bacteroidales
    ## Taxa_00234                           o__Cytophagales
    ## Taxa_00235                           o__Cytophagales
    ## Taxa_00236                           o__Cytophagales
    ## Taxa_00237                           o__Cytophagales
    ## Taxa_00238                           o__Cytophagales
    ## Taxa_00239                           o__Cytophagales
    ## Taxa_00240                           o__Cytophagales
    ## Taxa_00241                           o__Cytophagales
    ## Taxa_00242                           o__Cytophagales
    ## Taxa_00243                           o__Cytophagales
    ## Taxa_00244                           o__Cytophagales
    ## Taxa_00245                           o__Cytophagales
    ## Taxa_00246                           o__Cytophagales
    ## Taxa_00247                           o__Cytophagales
    ## Taxa_00248                           o__Cytophagales
    ## Taxa_00249                           o__Cytophagales
    ## Taxa_00250                           o__Cytophagales
    ## Taxa_00251                           o__Cytophagales
    ## Taxa_00252                           o__Cytophagales
    ## Taxa_00253                           o__Cytophagales
    ## Taxa_00254                           o__Cytophagales
    ## Taxa_00255                           o__Cytophagales
    ## Taxa_00256                           o__Cytophagales
    ## Taxa_00257                       o__Flavobacteriales
    ## Taxa_00258                       o__Flavobacteriales
    ## Taxa_00259                       o__Flavobacteriales
    ## Taxa_00260                       o__Flavobacteriales
    ## Taxa_00261                       o__Flavobacteriales
    ## Taxa_00262                       o__Flavobacteriales
    ## Taxa_00263                       o__Flavobacteriales
    ## Taxa_00264                       o__Flavobacteriales
    ## Taxa_00265                       o__Flavobacteriales
    ## Taxa_00266                       o__Flavobacteriales
    ## Taxa_00267                       o__Flavobacteriales
    ## Taxa_00268                       o__Flavobacteriales
    ## Taxa_00269                       o__Flavobacteriales
    ## Taxa_00270                     o__Sphingobacteriales
    ## Taxa_00271                     o__Sphingobacteriales
    ## Taxa_00272                     o__Sphingobacteriales
    ## Taxa_00273                     o__Sphingobacteriales
    ## Taxa_00274                     o__Sphingobacteriales
    ## Taxa_00275                     o__Sphingobacteriales
    ## Taxa_00276                     o__Sphingobacteriales
    ## Taxa_00277                     o__Sphingobacteriales
    ## Taxa_00278                     o__Sphingobacteriales
    ## Taxa_00279                     o__Sphingobacteriales
    ## Taxa_00280                     o__Sphingobacteriales
    ## Taxa_00281                     o__Sphingobacteriales
    ## Taxa_00282                     o__Sphingobacteriales
    ## Taxa_00283                     o__Sphingobacteriales
    ## Taxa_00284                     o__Sphingobacteriales
    ## Taxa_00285                     o__Sphingobacteriales
    ## Taxa_00286                     o__Sphingobacteriales
    ## Taxa_00287                     o__Sphingobacteriales
    ## Taxa_00288                     o__Sphingobacteriales
    ## Taxa_00289                     o__Sphingobacteriales
    ## Taxa_00290                     o__Sphingobacteriales
    ## Taxa_00291                     o__Sphingobacteriales
    ## Taxa_00292                     o__Sphingobacteriales
    ## Taxa_00293                     o__Sphingobacteriales
    ## Taxa_00294                     o__Sphingobacteriales
    ## Taxa_00295                     o__Sphingobacteriales
    ## Taxa_00296                     o__Sphingobacteriales
    ## Taxa_00297                     o__Sphingobacteriales
    ## Taxa_00298                     o__Sphingobacteriales
    ## Taxa_00299                     o__Sphingobacteriales
    ## Taxa_00300                     o__Sphingobacteriales
    ## Taxa_00301                     o__Sphingobacteriales
    ## Taxa_00302                     o__Sphingobacteriales
    ## Taxa_00303                     o__Sphingobacteriales
    ## Taxa_00304             o__BRC1_genera_incertae_sedis
    ## Taxa_00305            o__WPS-1_genera_incertae_sedis
    ## Taxa_00306            o__WPS-2_genera_incertae_sedis
    ## Taxa_00307              o__ZB3_genera_incertae_sedis
    ## Taxa_00308 o__Saccharibacteria_genera_incertae_sedis
    ## Taxa_00309                           o__Chlamydiales
    ## Taxa_00310                           o__Chlamydiales
    ## Taxa_00311                           o__Chlamydiales
    ## Taxa_00312                           o__Chlamydiales
    ## Taxa_00313                           o__Chlamydiales
    ## Taxa_00314                                      <NA>
    ## Taxa_00315                         o__Anaerolineales
    ## Taxa_00316                         o__Anaerolineales
    ## Taxa_00317                         o__Anaerolineales
    ## Taxa_00318                         o__Anaerolineales
    ## Taxa_00319                         o__Anaerolineales
    ## Taxa_00320                         o__Anaerolineales
    ## Taxa_00321                       o__Ardenticatenales
    ## Taxa_00322                          o__Caldilineales
    ## Taxa_00323                          o__Caldilineales
    ## Taxa_00324                          o__Caldilineales
    ## Taxa_00325                                      <NA>
    ## Taxa_00326                         o__Chloroflexales
    ## Taxa_00327                         o__Chloroflexales
    ## Taxa_00328                         o__Chloroflexales
    ## Taxa_00329                         o__Chloroflexales
    ## Taxa_00330                          o__Kallotenuales
    ## Taxa_00331                      o__Dehalococcoidales
    ## Taxa_00332                                      <NA>
    ## Taxa_00333                      o__Ktedonobacterales
    ## Taxa_00334                      o__Ktedonobacterales
    ## Taxa_00335                      o__Ktedonobacterales
    ## Taxa_00336                         o__Thermoflexales
    ## Taxa_00337                                      <NA>
    ## Taxa_00338                      o__Sphaerobacterales
    ## Taxa_00339                      o__Sphaerobacterales
    ## Taxa_00340                      o__Thermomicrobiales
    ## Taxa_00341                                      <NA>
    ## Taxa_00342                            o__Chloroplast
    ## Taxa_00343                            o__Chloroplast
    ## Taxa_00344                            o__Chloroplast
    ## Taxa_00345                            o__Chloroplast
    ## Taxa_00346                            o__Chloroplast
    ## Taxa_00347                            o__Chloroplast
    ## Taxa_00348                               o__Family I
    ## Taxa_00349                                      <NA>
    ## Taxa_00350                                      <NA>
    ## Taxa_00351                       o__Elusimicrobiales
    ## Taxa_00352               o__Candidatus Endomicrobium
    ## Taxa_00353                                      <NA>
    ## Taxa_00354                                      <NA>
    ## Taxa_00355                             o__Bacillales
    ## Taxa_00356                             o__Bacillales
    ## Taxa_00357                             o__Bacillales
    ## Taxa_00358                             o__Bacillales
    ## Taxa_00359                             o__Bacillales
    ## Taxa_00360                             o__Bacillales
    ## Taxa_00361                             o__Bacillales
    ## Taxa_00362                             o__Bacillales
    ## Taxa_00363                             o__Bacillales
    ## Taxa_00364                             o__Bacillales
    ## Taxa_00365                             o__Bacillales
    ## Taxa_00366                             o__Bacillales
    ## Taxa_00367                             o__Bacillales
    ## Taxa_00368                             o__Bacillales
    ## Taxa_00369                             o__Bacillales
    ## Taxa_00370                             o__Bacillales
    ## Taxa_00371                             o__Bacillales
    ## Taxa_00372                             o__Bacillales
    ## Taxa_00373                             o__Bacillales
    ## Taxa_00374                             o__Bacillales
    ## Taxa_00375                             o__Bacillales
    ## Taxa_00376                             o__Bacillales
    ## Taxa_00377                             o__Bacillales
    ## Taxa_00378                             o__Bacillales
    ## Taxa_00379                             o__Bacillales
    ## Taxa_00380                             o__Bacillales
    ## Taxa_00381                             o__Bacillales
    ## Taxa_00382                             o__Bacillales
    ## Taxa_00383                             o__Bacillales
    ## Taxa_00384                             o__Bacillales
    ## Taxa_00385                             o__Bacillales
    ## Taxa_00386                             o__Bacillales
    ## Taxa_00387                             o__Bacillales
    ## Taxa_00388                             o__Bacillales
    ## Taxa_00389                             o__Bacillales
    ## Taxa_00390                             o__Bacillales
    ## Taxa_00391                             o__Bacillales
    ## Taxa_00392                             o__Bacillales
    ## Taxa_00393                             o__Bacillales
    ## Taxa_00394                             o__Bacillales
    ## Taxa_00395                             o__Bacillales
    ## Taxa_00396                             o__Bacillales
    ## Taxa_00397                             o__Bacillales
    ## Taxa_00398                             o__Bacillales
    ## Taxa_00399                             o__Bacillales
    ## Taxa_00400                             o__Bacillales
    ## Taxa_00401                             o__Bacillales
    ## Taxa_00402                             o__Bacillales
    ## Taxa_00403                             o__Bacillales
    ## Taxa_00404                             o__Bacillales
    ## Taxa_00405                             o__Bacillales
    ## Taxa_00406                             o__Bacillales
    ## Taxa_00407                             o__Bacillales
    ## Taxa_00408                             o__Bacillales
    ## Taxa_00409                             o__Bacillales
    ## Taxa_00410                             o__Bacillales
    ## Taxa_00411                             o__Bacillales
    ## Taxa_00412                             o__Bacillales
    ## Taxa_00413                             o__Bacillales
    ## Taxa_00414                             o__Bacillales
    ## Taxa_00415                             o__Bacillales
    ## Taxa_00416                             o__Bacillales
    ## Taxa_00417                             o__Bacillales
    ## Taxa_00418                        o__Lactobacillales
    ## Taxa_00419                        o__Lactobacillales
    ## Taxa_00420                        o__Lactobacillales
    ## Taxa_00421                        o__Lactobacillales
    ## Taxa_00422                        o__Lactobacillales
    ## Taxa_00423                        o__Lactobacillales
    ## Taxa_00424                                      <NA>
    ## Taxa_00425                          o__Clostridiales
    ## Taxa_00426                          o__Clostridiales
    ## Taxa_00427                          o__Clostridiales
    ## Taxa_00428                          o__Clostridiales
    ## Taxa_00429                          o__Clostridiales
    ## Taxa_00430                          o__Clostridiales
    ## Taxa_00431                          o__Clostridiales
    ## Taxa_00432                          o__Clostridiales
    ## Taxa_00433                          o__Clostridiales
    ## Taxa_00434                          o__Clostridiales
    ## Taxa_00435                          o__Clostridiales
    ## Taxa_00436                          o__Clostridiales
    ## Taxa_00437                          o__Clostridiales
    ## Taxa_00438                          o__Clostridiales
    ## Taxa_00439                          o__Clostridiales
    ## Taxa_00440                          o__Clostridiales
    ## Taxa_00441                          o__Clostridiales
    ## Taxa_00442                          o__Clostridiales
    ## Taxa_00443                          o__Clostridiales
    ## Taxa_00444                          o__Clostridiales
    ## Taxa_00445                          o__Clostridiales
    ## Taxa_00446                          o__Clostridiales
    ## Taxa_00447                          o__Clostridiales
    ## Taxa_00448                          o__Clostridiales
    ## Taxa_00449                          o__Clostridiales
    ## Taxa_00450                          o__Clostridiales
    ## Taxa_00451                          o__Clostridiales
    ## Taxa_00452                          o__Clostridiales
    ## Taxa_00453                          o__Clostridiales
    ## Taxa_00454                          o__Clostridiales
    ## Taxa_00455                          o__Clostridiales
    ## Taxa_00456                          o__Clostridiales
    ## Taxa_00457                          o__Clostridiales
    ## Taxa_00458                          o__Clostridiales
    ## Taxa_00459                          o__Clostridiales
    ## Taxa_00460                 o__Thermoanaerobacterales
    ## Taxa_00461                     o__Erysipelotrichales
    ## Taxa_00462                        o__Selenomonadales
    ## Taxa_00463                        o__Selenomonadales
    ## Taxa_00464                        o__Selenomonadales
    ## Taxa_00465                        o__Selenomonadales
    ## Taxa_00466                        o__Selenomonadales
    ## Taxa_00467                       o__Gemmatimonadales
    ## Taxa_00468               o__Candidatus Hydrogenedens
    ## Taxa_00469  o__Latescibacteria_genera_incertae_sedis
    ## Taxa_00470                                      <NA>
    ## Taxa_00471                          o__Victivallales
    ## Taxa_00472   o__Microgenomates_genera_incertae_sedis
    ## Taxa_00473                          o__Nitrospirales
    ## Taxa_00474     o__Omnitrophica_genera_incertae_sedis
    ## Taxa_00475    o__Parcubacteria_genera_incertae_sedis
    ## Taxa_00476                                      <NA>
    ## Taxa_00477                                      <NA>
    ## Taxa_00478                        o__Phycisphaerales
    ## Taxa_00479                       o__Tepidisphaerales
    ## Taxa_00480                                      <NA>
    ## Taxa_00481                 o__Candidatus Brocadiales
    ## Taxa_00482                 o__Candidatus Brocadiales
    ## Taxa_00483                       o__Planctomycetales
    ## Taxa_00484                       o__Planctomycetales
    ## Taxa_00485                       o__Planctomycetales
    ## Taxa_00486                       o__Planctomycetales
    ## Taxa_00487                       o__Planctomycetales
    ## Taxa_00488                       o__Planctomycetales
    ## Taxa_00489                       o__Planctomycetales
    ## Taxa_00490                       o__Planctomycetales
    ## Taxa_00491                       o__Planctomycetales
    ## Taxa_00492                       o__Planctomycetales
    ## Taxa_00493                       o__Planctomycetales
    ## Taxa_00494                       o__Planctomycetales
    ## Taxa_00495                       o__Planctomycetales
    ## Taxa_00496                       o__Planctomycetales
    ## Taxa_00497                       o__Planctomycetales
    ## Taxa_00498                                      <NA>
    ## Taxa_00499                                      <NA>
    ## Taxa_00500     o__Alphaproteobacteria_incertae_sedis
    ## Taxa_00501     o__Alphaproteobacteria_incertae_sedis
    ## Taxa_00502     o__Alphaproteobacteria_incertae_sedis
    ## Taxa_00503     o__Alphaproteobacteria_incertae_sedis
    ## Taxa_00504                        o__Caulobacterales
    ## Taxa_00505                        o__Caulobacterales
    ## Taxa_00506                        o__Caulobacterales
    ## Taxa_00507                        o__Caulobacterales
    ## Taxa_00508                        o__Caulobacterales
    ## Taxa_00509                        o__Caulobacterales
    ## Taxa_00510                        o__Caulobacterales
    ## Taxa_00511                        o__Caulobacterales
    ## Taxa_00512                        o__Caulobacterales
    ## Taxa_00513                            o__Eilatimonas
    ## Taxa_00514                        o__Parvularculales
    ## Taxa_00515                            o__Rhizobiales
    ## Taxa_00516                            o__Rhizobiales
    ## Taxa_00517                            o__Rhizobiales
    ## Taxa_00518                            o__Rhizobiales
    ## Taxa_00519                            o__Rhizobiales
    ## Taxa_00520                            o__Rhizobiales
    ## Taxa_00521                            o__Rhizobiales
    ## Taxa_00522                            o__Rhizobiales
    ## Taxa_00523                            o__Rhizobiales
    ## Taxa_00524                            o__Rhizobiales
    ## Taxa_00525                            o__Rhizobiales
    ## Taxa_00526                            o__Rhizobiales
    ## Taxa_00527                            o__Rhizobiales
    ## Taxa_00528                            o__Rhizobiales
    ## Taxa_00529                            o__Rhizobiales
    ## Taxa_00530                            o__Rhizobiales
    ## Taxa_00531                            o__Rhizobiales
    ## Taxa_00532                            o__Rhizobiales
    ## Taxa_00533                            o__Rhizobiales
    ## Taxa_00534                            o__Rhizobiales
    ## Taxa_00535                            o__Rhizobiales
    ## Taxa_00536                            o__Rhizobiales
    ## Taxa_00537                            o__Rhizobiales
    ## Taxa_00538                            o__Rhizobiales
    ## Taxa_00539                            o__Rhizobiales
    ## Taxa_00540                            o__Rhizobiales
    ## Taxa_00541                            o__Rhizobiales
    ## Taxa_00542                            o__Rhizobiales
    ## Taxa_00543                            o__Rhizobiales
    ## Taxa_00544                            o__Rhizobiales
    ## Taxa_00545                            o__Rhizobiales
    ## Taxa_00546                            o__Rhizobiales
    ## Taxa_00547                            o__Rhizobiales
    ## Taxa_00548                            o__Rhizobiales
    ## Taxa_00549                            o__Rhizobiales
    ## Taxa_00550                            o__Rhizobiales
    ## Taxa_00551                            o__Rhizobiales
    ## Taxa_00552                            o__Rhizobiales
    ## Taxa_00553                            o__Rhizobiales
    ## Taxa_00554                            o__Rhizobiales
    ## Taxa_00555                            o__Rhizobiales
    ## Taxa_00556                            o__Rhizobiales
    ## Taxa_00557                            o__Rhizobiales
    ## Taxa_00558                            o__Rhizobiales
    ## Taxa_00559                            o__Rhizobiales
    ## Taxa_00560                            o__Rhizobiales
    ## Taxa_00561                            o__Rhizobiales
    ## Taxa_00562                            o__Rhizobiales
    ## Taxa_00563                            o__Rhizobiales
    ## Taxa_00564                            o__Rhizobiales
    ## Taxa_00565                            o__Rhizobiales
    ## Taxa_00566                            o__Rhizobiales
    ## Taxa_00567                            o__Rhizobiales
    ## Taxa_00568                            o__Rhizobiales
    ## Taxa_00569                            o__Rhizobiales
    ## Taxa_00570                            o__Rhizobiales
    ## Taxa_00571                            o__Rhizobiales
    ## Taxa_00572                            o__Rhizobiales
    ## Taxa_00573                            o__Rhizobiales
    ## Taxa_00574                            o__Rhizobiales
    ## Taxa_00575                            o__Rhizobiales
    ## Taxa_00576                            o__Rhizobiales
    ## Taxa_00577                            o__Rhizobiales
    ## Taxa_00578                            o__Rhizobiales
    ## Taxa_00579                            o__Rhizobiales
    ## Taxa_00580                            o__Rhizobiales
    ## Taxa_00581                            o__Rhizobiales
    ## Taxa_00582                            o__Rhizobiales
    ## Taxa_00583                            o__Rhizobiales
    ## Taxa_00584                            o__Rhizobiales
    ## Taxa_00585                            o__Rhizobiales
    ## Taxa_00586                            o__Rhizobiales
    ## Taxa_00587                            o__Rhizobiales
    ## Taxa_00588                            o__Rhizobiales
    ## Taxa_00589                            o__Rhizobiales
    ## Taxa_00590                            o__Rhizobiales
    ## Taxa_00591                            o__Rhizobiales
    ## Taxa_00592                            o__Rhizobiales
    ## Taxa_00593                            o__Rhizobiales
    ## Taxa_00594                            o__Rhizobiales
    ## Taxa_00595                            o__Rhizobiales
    ## Taxa_00596                            o__Rhizobiales
    ## Taxa_00597                            o__Rhizobiales
    ## Taxa_00598                            o__Rhizobiales
    ## Taxa_00599                        o__Rhodobacterales
    ## Taxa_00600                        o__Rhodobacterales
    ## Taxa_00601                        o__Rhodobacterales
    ## Taxa_00602                        o__Rhodobacterales
    ## Taxa_00603                        o__Rhodobacterales
    ## Taxa_00604                        o__Rhodobacterales
    ## Taxa_00605                        o__Rhodobacterales
    ## Taxa_00606                        o__Rhodobacterales
    ## Taxa_00607                       o__Rhodospirillales
    ## Taxa_00608                       o__Rhodospirillales
    ## Taxa_00609                       o__Rhodospirillales
    ## Taxa_00610                       o__Rhodospirillales
    ## Taxa_00611                       o__Rhodospirillales
    ## Taxa_00612                       o__Rhodospirillales
    ## Taxa_00613                       o__Rhodospirillales
    ## Taxa_00614                       o__Rhodospirillales
    ## Taxa_00615                       o__Rhodospirillales
    ## Taxa_00616                       o__Rhodospirillales
    ## Taxa_00617                       o__Rhodospirillales
    ## Taxa_00618                       o__Rhodospirillales
    ## Taxa_00619                       o__Rhodospirillales
    ## Taxa_00620                       o__Rhodospirillales
    ## Taxa_00621                       o__Rhodospirillales
    ## Taxa_00622                       o__Rhodospirillales
    ## Taxa_00623                       o__Rhodospirillales
    ## Taxa_00624                       o__Rhodospirillales
    ## Taxa_00625                       o__Rhodospirillales
    ## Taxa_00626                       o__Rhodospirillales
    ## Taxa_00627                       o__Rhodospirillales
    ## Taxa_00628                       o__Rhodospirillales
    ## Taxa_00629                       o__Rhodospirillales
    ## Taxa_00630                       o__Rhodospirillales
    ## Taxa_00631                       o__Rhodospirillales
    ## Taxa_00632                       o__Rhodospirillales
    ## Taxa_00633                       o__Rhodospirillales
    ## Taxa_00634                       o__Rhodospirillales
    ## Taxa_00635                       o__Rhodospirillales
    ## Taxa_00636                       o__Rhodospirillales
    ## Taxa_00637                       o__Rhodospirillales
    ## Taxa_00638                       o__Rhodospirillales
    ## Taxa_00639                       o__Rhodospirillales
    ## Taxa_00640                       o__Rhodospirillales
    ## Taxa_00641                       o__Rhodospirillales
    ## Taxa_00642                       o__Rhodospirillales
    ## Taxa_00643                       o__Rhodospirillales
    ## Taxa_00644                       o__Rhodospirillales
    ## Taxa_00645                       o__Rhodospirillales
    ## Taxa_00646                       o__Rhodospirillales
    ## Taxa_00647                       o__Rhodospirillales
    ## Taxa_00648                       o__Rhodospirillales
    ## Taxa_00649                       o__Rhodospirillales
    ## Taxa_00650                       o__Rhodospirillales
    ## Taxa_00651                       o__Rhodospirillales
    ## Taxa_00652                       o__Rhodospirillales
    ## Taxa_00653                          o__Rickettsiales
    ## Taxa_00654                          o__Rickettsiales
    ## Taxa_00655                          o__Rickettsiales
    ## Taxa_00656                          o__Rickettsiales
    ## Taxa_00657                          o__Rickettsiales
    ## Taxa_00658                         o__Sneathiellales
    ## Taxa_00659                         o__Sneathiellales
    ## Taxa_00660                       o__Sphingomonadales
    ## Taxa_00661                       o__Sphingomonadales
    ## Taxa_00662                       o__Sphingomonadales
    ## Taxa_00663                       o__Sphingomonadales
    ## Taxa_00664                       o__Sphingomonadales
    ## Taxa_00665                       o__Sphingomonadales
    ## Taxa_00666                       o__Sphingomonadales
    ## Taxa_00667                       o__Sphingomonadales
    ## Taxa_00668                       o__Sphingomonadales
    ## Taxa_00669                       o__Sphingomonadales
    ## Taxa_00670                       o__Sphingomonadales
    ## Taxa_00671                       o__Sphingomonadales
    ## Taxa_00672                       o__Sphingomonadales
    ## Taxa_00673                       o__Sphingomonadales
    ## Taxa_00674                       o__Sphingomonadales
    ## Taxa_00675                       o__Sphingomonadales
    ## Taxa_00676                       o__Sphingomonadales
    ## Taxa_00677                       o__Sphingomonadales
    ## Taxa_00678                       o__Sphingomonadales
    ## Taxa_00679                       o__Sphingomonadales
    ## Taxa_00680                       o__Sphingomonadales
    ## Taxa_00681                       o__Sphingomonadales
    ## Taxa_00682                                      <NA>
    ## Taxa_00683                        o__Burkholderiales
    ## Taxa_00684                        o__Burkholderiales
    ## Taxa_00685                        o__Burkholderiales
    ## Taxa_00686                        o__Burkholderiales
    ## Taxa_00687                        o__Burkholderiales
    ## Taxa_00688                        o__Burkholderiales
    ## Taxa_00689                        o__Burkholderiales
    ## Taxa_00690                        o__Burkholderiales
    ## Taxa_00691                        o__Burkholderiales
    ## Taxa_00692                        o__Burkholderiales
    ## Taxa_00693                        o__Burkholderiales
    ## Taxa_00694                        o__Burkholderiales
    ## Taxa_00695                        o__Burkholderiales
    ## Taxa_00696                        o__Burkholderiales
    ## Taxa_00697                        o__Burkholderiales
    ## Taxa_00698                        o__Burkholderiales
    ## Taxa_00699                        o__Burkholderiales
    ## Taxa_00700                        o__Burkholderiales
    ## Taxa_00701                        o__Burkholderiales
    ## Taxa_00702                        o__Burkholderiales
    ## Taxa_00703                        o__Burkholderiales
    ## Taxa_00704                        o__Burkholderiales
    ## Taxa_00705                        o__Burkholderiales
    ## Taxa_00706                        o__Burkholderiales
    ## Taxa_00707                        o__Burkholderiales
    ## Taxa_00708                        o__Burkholderiales
    ## Taxa_00709                        o__Burkholderiales
    ## Taxa_00710                        o__Burkholderiales
    ## Taxa_00711                        o__Burkholderiales
    ## Taxa_00712                        o__Burkholderiales
    ## Taxa_00713                        o__Burkholderiales
    ## Taxa_00714                        o__Burkholderiales
    ## Taxa_00715                        o__Burkholderiales
    ## Taxa_00716                        o__Burkholderiales
    ## Taxa_00717                        o__Burkholderiales
    ## Taxa_00718                        o__Burkholderiales
    ## Taxa_00719                        o__Burkholderiales
    ## Taxa_00720                        o__Burkholderiales
    ## Taxa_00721                        o__Burkholderiales
    ## Taxa_00722                        o__Burkholderiales
    ## Taxa_00723                        o__Burkholderiales
    ## Taxa_00724                        o__Burkholderiales
    ## Taxa_00725                        o__Burkholderiales
    ## Taxa_00726                        o__Burkholderiales
    ## Taxa_00727                        o__Burkholderiales
    ## Taxa_00728                        o__Burkholderiales
    ## Taxa_00729                        o__Burkholderiales
    ## Taxa_00730                        o__Burkholderiales
    ## Taxa_00731                        o__Burkholderiales
    ## Taxa_00732                        o__Burkholderiales
    ## Taxa_00733                        o__Burkholderiales
    ## Taxa_00734                        o__Burkholderiales
    ## Taxa_00735                        o__Burkholderiales
    ## Taxa_00736                        o__Burkholderiales
    ## Taxa_00737                        o__Burkholderiales
    ## Taxa_00738                        o__Burkholderiales
    ## Taxa_00739                        o__Burkholderiales
    ## Taxa_00740                        o__Burkholderiales
    ## Taxa_00741                        o__Burkholderiales
    ## Taxa_00742                        o__Burkholderiales
    ## Taxa_00743                        o__Burkholderiales
    ## Taxa_00744                        o__Burkholderiales
    ## Taxa_00745                        o__Burkholderiales
    ## Taxa_00746                        o__Burkholderiales
    ## Taxa_00747                        o__Burkholderiales
    ## Taxa_00748                        o__Burkholderiales
    ## Taxa_00749                        o__Burkholderiales
    ## Taxa_00750                       o__Ferritrophicales
    ## Taxa_00751                             o__Ferrovales
    ## Taxa_00752                         o__Gallionellales
    ## Taxa_00753                         o__Gallionellales
    ## Taxa_00754                         o__Gallionellales
    ## Taxa_00755                        o__Methylophilales
    ## Taxa_00756                        o__Methylophilales
    ## Taxa_00757                           o__Neisseriales
    ## Taxa_00758                           o__Neisseriales
    ## Taxa_00759                           o__Neisseriales
    ## Taxa_00760                           o__Neisseriales
    ## Taxa_00761                           o__Neisseriales
    ## Taxa_00762                           o__Neisseriales
    ## Taxa_00763                           o__Neisseriales
    ## Taxa_00764                       o__Nitrosomonadales
    ## Taxa_00765                       o__Nitrosomonadales
    ## Taxa_00766                       o__Procabacteriales
    ## Taxa_00767                          o__Rhodocyclales
    ## Taxa_00768                          o__Rhodocyclales
    ## Taxa_00769                          o__Rhodocyclales
    ## Taxa_00770                          o__Rhodocyclales
    ## Taxa_00771                          o__Rhodocyclales
    ## Taxa_00772                          o__Rhodocyclales
    ## Taxa_00773                          o__Rhodocyclales
    ## Taxa_00774                          o__Rhodocyclales
    ## Taxa_00775                          o__Rhodocyclales
    ## Taxa_00776                          o__Rhodocyclales
    ## Taxa_00777                        o__Sulfuricellales
    ## Taxa_00778                        o__Sulfuricellales
    ## Taxa_00779                                      <NA>
    ## Taxa_00780                      o__Bdellovibrionales
    ## Taxa_00781                      o__Bdellovibrionales
    ## Taxa_00782                      o__Bdellovibrionales
    ## Taxa_00783                      o__Bdellovibrionales
    ## Taxa_00784                      o__Bdellovibrionales
    ## Taxa_00785                      o__Bdellovibrionales
    ## Taxa_00786                      o__Bdellovibrionales
    ## Taxa_00787                      o__Bdellovibrionales
    ## Taxa_00788     o__Deltaproteobacteria_incertae_sedis
    ## Taxa_00789     o__Deltaproteobacteria_incertae_sedis
    ## Taxa_00790     o__Deltaproteobacteria_incertae_sedis
    ## Taxa_00791                      o__Desulfobacterales
    ## Taxa_00792                      o__Desulfobacterales
    ## Taxa_00793                      o__Desulfobacterales
    ## Taxa_00794                     o__Desulfovibrionales
    ## Taxa_00795                     o__Desulfovibrionales
    ## Taxa_00796                     o__Desulfuromonadales
    ## Taxa_00797                     o__Desulfuromonadales
    ## Taxa_00798                     o__Desulfuromonadales
    ## Taxa_00799                     o__Desulfuromonadales
    ## Taxa_00800                           o__Myxococcales
    ## Taxa_00801                           o__Myxococcales
    ## Taxa_00802                           o__Myxococcales
    ## Taxa_00803                           o__Myxococcales
    ## Taxa_00804                           o__Myxococcales
    ## Taxa_00805                           o__Myxococcales
    ## Taxa_00806                           o__Myxococcales
    ## Taxa_00807                           o__Myxococcales
    ## Taxa_00808                           o__Myxococcales
    ## Taxa_00809                           o__Myxococcales
    ## Taxa_00810                           o__Myxococcales
    ## Taxa_00811                           o__Myxococcales
    ## Taxa_00812                           o__Myxococcales
    ## Taxa_00813                           o__Myxococcales
    ## Taxa_00814                           o__Myxococcales
    ## Taxa_00815                           o__Myxococcales
    ## Taxa_00816                           o__Myxococcales
    ## Taxa_00817                           o__Myxococcales
    ## Taxa_00818                           o__Myxococcales
    ## Taxa_00819                           o__Myxococcales
    ## Taxa_00820                           o__Myxococcales
    ## Taxa_00821                           o__Myxococcales
    ## Taxa_00822                           o__Myxococcales
    ## Taxa_00823                           o__Myxococcales
    ## Taxa_00824                           o__Myxococcales
    ## Taxa_00825                           o__Myxococcales
    ## Taxa_00826                           o__Myxococcales
    ## Taxa_00827                    o__Syntrophobacterales
    ## Taxa_00828                    o__Syntrophobacterales
    ## Taxa_00829                    o__Syntrophobacterales
    ## Taxa_00830                                      <NA>
    ## Taxa_00831                            o__Nautiliales
    ## Taxa_00832                                      <NA>
    ## Taxa_00833                           o__Chromatiales
    ## Taxa_00834                           o__Chromatiales
    ## Taxa_00835                           o__Chromatiales
    ## Taxa_00836                           o__Chromatiales
    ## Taxa_00837                           o__Chromatiales
    ## Taxa_00838                           o__Chromatiales
    ## Taxa_00839                           o__Chromatiales
    ## Taxa_00840                           o__Chromatiales
    ## Taxa_00841                           o__Chromatiales
    ## Taxa_00842                           o__Chromatiales
    ## Taxa_00843                           o__Chromatiales
    ## Taxa_00844                           o__Chromatiales
    ## Taxa_00845                           o__Chromatiales
    ## Taxa_00846                           o__Chromatiales
    ## Taxa_00847                           o__Chromatiales
    ## Taxa_00848                      o__Enterobacteriales
    ## Taxa_00849                      o__Enterobacteriales
    ## Taxa_00850                      o__Enterobacteriales
    ## Taxa_00851                      o__Enterobacteriales
    ## Taxa_00852                      o__Enterobacteriales
    ## Taxa_00853     o__Gammaproteobacteria_incertae_sedis
    ## Taxa_00854     o__Gammaproteobacteria_incertae_sedis
    ## Taxa_00855     o__Gammaproteobacteria_incertae_sedis
    ## Taxa_00856     o__Gammaproteobacteria_incertae_sedis
    ## Taxa_00857                          o__Legionellales
    ## Taxa_00858                          o__Legionellales
    ## Taxa_00859                          o__Legionellales
    ## Taxa_00860                          o__Legionellales
    ## Taxa_00861                          o__Legionellales
    ## Taxa_00862                          o__Legionellales
    ## Taxa_00863                        o__Methylococcales
    ## Taxa_00864                        o__Methylococcales
    ## Taxa_00865                        o__Methylococcales
    ## Taxa_00866                        o__Methylococcales
    ## Taxa_00867                        o__Methylococcales
    ## Taxa_00868                      o__Oceanospirillales
    ## Taxa_00869                      o__Oceanospirillales
    ## Taxa_00870                      o__Oceanospirillales
    ## Taxa_00871                      o__Oceanospirillales
    ## Taxa_00872                         o__Pasteurellales
    ## Taxa_00873                        o__Pseudomonadales
    ## Taxa_00874                        o__Pseudomonadales
    ## Taxa_00875                        o__Pseudomonadales
    ## Taxa_00876                        o__Pseudomonadales
    ## Taxa_00877                        o__Pseudomonadales
    ## Taxa_00878                        o__Pseudomonadales
    ## Taxa_00879                        o__Pseudomonadales
    ## Taxa_00880                        o__Pseudomonadales
    ## Taxa_00881                          o__Thiotrichales
    ## Taxa_00882                          o__Thiotrichales
    ## Taxa_00883                        o__Xanthomonadales
    ## Taxa_00884                        o__Xanthomonadales
    ## Taxa_00885                        o__Xanthomonadales
    ## Taxa_00886                        o__Xanthomonadales
    ## Taxa_00887                        o__Xanthomonadales
    ## Taxa_00888                        o__Xanthomonadales
    ## Taxa_00889                        o__Xanthomonadales
    ## Taxa_00890                        o__Xanthomonadales
    ## Taxa_00891                        o__Xanthomonadales
    ## Taxa_00892                        o__Xanthomonadales
    ## Taxa_00893                        o__Xanthomonadales
    ## Taxa_00894                        o__Xanthomonadales
    ## Taxa_00895                        o__Xanthomonadales
    ## Taxa_00896                        o__Xanthomonadales
    ## Taxa_00897                        o__Xanthomonadales
    ## Taxa_00898                        o__Xanthomonadales
    ## Taxa_00899                        o__Xanthomonadales
    ## Taxa_00900                        o__Xanthomonadales
    ## Taxa_00901                        o__Xanthomonadales
    ## Taxa_00902                        o__Xanthomonadales
    ## Taxa_00903                        o__Xanthomonadales
    ## Taxa_00904                        o__Xanthomonadales
    ## Taxa_00905                        o__Xanthomonadales
    ## Taxa_00906                        o__Xanthomonadales
    ## Taxa_00907                        o__Xanthomonadales
    ## Taxa_00908                        o__Xanthomonadales
    ## Taxa_00909                        o__Xanthomonadales
    ## Taxa_00910                        o__Xanthomonadales
    ## Taxa_00911                        o__Xanthomonadales
    ## Taxa_00912                        o__Xanthomonadales
    ## Taxa_00913                        o__Xanthomonadales
    ## Taxa_00914                          o__Oligoflexales
    ## Taxa_00915                         o__Spirochaetales
    ## Taxa_00916                         o__Spirochaetales
    ## Taxa_00917                         o__Spirochaetales
    ## Taxa_00918                         o__Spirochaetales
    ## Taxa_00919                         o__Spirochaetales
    ## Taxa_00920                         o__Spirochaetales
    ## Taxa_00921                         o__Spirochaetales
    ## Taxa_00922                         o__Spirochaetales
    ## Taxa_00923              o__SR1_genera_incertae_sedis
    ## Taxa_00924                      o__Entomoplasmatales
    ## Taxa_00925                        o__Haloplasmatales
    ## Taxa_00926               o__Thermodesulfobacteriales
    ## Taxa_00927                                      <NA>
    ## Taxa_00928                                      <NA>
    ## Taxa_00929                             o__Opitutales
    ## Taxa_00930                             o__Opitutales
    ## Taxa_00931                             o__Opitutales
    ## Taxa_00932                        o__Puniceicoccales
    ## Taxa_00933                        o__Puniceicoccales
    ## Taxa_00934                        o__Puniceicoccales
    ## Taxa_00935                                      <NA>
    ## Taxa_00936   o__Spartobacteria_genera_incertae_sedis
    ## Taxa_00937                         o__Terrimicrobium
    ## Taxa_00938                                      <NA>
    ## Taxa_00939                            o__Limisphaera
    ## Taxa_00940     o__Subdivision3_genera_incertae_sedis
    ## Taxa_00941                     o__Verrucomicrobiales
    ## Taxa_00942                     o__Verrucomicrobiales
    ## Taxa_00943                     o__Verrucomicrobiales
    ## Taxa_00944                     o__Verrucomicrobiales
    ## Taxa_00945                     o__Verrucomicrobiales
    ## Taxa_00946                     o__Verrucomicrobiales
    ## Taxa_00947                     o__Verrucomicrobiales
    ## Taxa_00948                     o__Verrucomicrobiales
    ## Taxa_00949                                      <NA>
    ##                                               Family
    ## Taxa_00000                                      <NA>
    ## Taxa_00001                         f__Pyrodictiaceae
    ## Taxa_00002                                      <NA>
    ## Taxa_00003                                      <NA>
    ## Taxa_00004                                      <NA>
    ## Taxa_00005                                   f__Gp10
    ## Taxa_00006                                   f__Gp11
    ## Taxa_00007                                   f__Gp12
    ## Taxa_00008                                   f__Gp13
    ## Taxa_00009                                   f__Gp15
    ## Taxa_00010                                   f__Gp16
    ## Taxa_00011                                   f__Gp17
    ## Taxa_00012                                   f__Gp18
    ## Taxa_00013                                   f__Gp19
    ## Taxa_00014                             f__Acidicapsa
    ## Taxa_00015                              f__Acidipila
    ## Taxa_00016                         f__Acidobacterium
    ## Taxa_00017                              f__Bryocella
    ## Taxa_00018                  f__Candidatus Koribacter
    ## Taxa_00019                           f__Edaphobacter
    ## Taxa_00020                                    f__Gp1
    ## Taxa_00021                           f__Granulicella
    ## Taxa_00022                          f__Telmatobacter
    ## Taxa_00023                            f__Terriglobus
    ## Taxa_00024                                   f__Gp20
    ## Taxa_00025                                   f__Gp22
    ## Taxa_00026                                   f__Gp25
    ## Taxa_00027                                    f__Gp2
    ## Taxa_00028                                      <NA>
    ## Taxa_00029                             f__Bryobacter
    ## Taxa_00030                  f__Candidatus Solibacter
    ## Taxa_00031                                    f__Gp3
    ## Taxa_00032                          f__Paludibaculum
    ## Taxa_00033                                      <NA>
    ## Taxa_00034                            f__Aridibacter
    ## Taxa_00035                          f__Blastocatella
    ## Taxa_00036                                    f__Gp4
    ## Taxa_00037                            f__Pyrinomonas
    ## Taxa_00038                                    f__Gp5
    ## Taxa_00039                                    f__Gp6
    ## Taxa_00040                                    f__Gp7
    ## Taxa_00041                          f__Holophagaceae
    ## Taxa_00042                          f__Holophagaceae
    ## Taxa_00043                          f__Holophagaceae
    ## Taxa_00044                                      <NA>
    ## Taxa_00045                                      <NA>
    ## Taxa_00046                                      <NA>
    ## Taxa_00047                      f__Acidimicrobiaceae
    ## Taxa_00048                      f__Acidimicrobiaceae
    ## Taxa_00049                      f__Acidimicrobiaceae
    ## Taxa_00050                      f__Acidimicrobiaceae
    ## Taxa_00051        f__Acidimicrobineae_incertae_sedis
    ## Taxa_00052                              f__Iamiaceae
    ## Taxa_00053                              f__Iamiaceae
    ## Taxa_00054                              f__Iamiaceae
    ## Taxa_00055                                      <NA>
    ## Taxa_00056                        f__Acidothermaceae
    ## Taxa_00057                        f__Actinospicaceae
    ## Taxa_00058                       f__Beutenbergiaceae
    ## Taxa_00059                       f__Beutenbergiaceae
    ## Taxa_00060                      f__Catenulisporaceae
    ## Taxa_00061                      f__Cellulomonadaceae
    ## Taxa_00062                      f__Cellulomonadaceae
    ## Taxa_00063                     f__Corynebacteriaceae
    ## Taxa_00064       f__Corynebacterineae_incertae_sedis
    ## Taxa_00065                    f__Cryptosporangiaceae
    ## Taxa_00066                    f__Cryptosporangiaceae
    ## Taxa_00067                    f__Cryptosporangiaceae
    ## Taxa_00068                    f__Cryptosporangiaceae
    ## Taxa_00069                          f__Demequinaceae
    ## Taxa_00070                       f__Dermabacteraceae
    ## Taxa_00071                       f__Dermabacteraceae
    ## Taxa_00072                         f__Dermacoccaceae
    ## Taxa_00073                       f__Dermatophilaceae
    ## Taxa_00074                       f__Dermatophilaceae
    ## Taxa_00075                    f__Geodermatophilaceae
    ## Taxa_00076                    f__Geodermatophilaceae
    ## Taxa_00077                    f__Geodermatophilaceae
    ## Taxa_00078                     f__Intrasporangiaceae
    ## Taxa_00079                     f__Intrasporangiaceae
    ## Taxa_00080                     f__Intrasporangiaceae
    ## Taxa_00081                     f__Intrasporangiaceae
    ## Taxa_00082                     f__Intrasporangiaceae
    ## Taxa_00083                        f__Kineosporiaceae
    ## Taxa_00084                        f__Kineosporiaceae
    ## Taxa_00085                        f__Kineosporiaceae
    ## Taxa_00086                        f__Kineosporiaceae
    ## Taxa_00087                      f__Microbacteriaceae
    ## Taxa_00088                      f__Microbacteriaceae
    ## Taxa_00089                      f__Microbacteriaceae
    ## Taxa_00090                      f__Microbacteriaceae
    ## Taxa_00091                      f__Microbacteriaceae
    ## Taxa_00092                      f__Microbacteriaceae
    ## Taxa_00093                      f__Microbacteriaceae
    ## Taxa_00094                      f__Microbacteriaceae
    ## Taxa_00095                      f__Microbacteriaceae
    ## Taxa_00096                      f__Microbacteriaceae
    ## Taxa_00097                      f__Microbacteriaceae
    ## Taxa_00098                      f__Microbacteriaceae
    ## Taxa_00099                      f__Microbacteriaceae
    ## Taxa_00100                      f__Microbacteriaceae
    ## Taxa_00101                      f__Microbacteriaceae
    ## Taxa_00102                      f__Microbacteriaceae
    ## Taxa_00103                      f__Microbacteriaceae
    ## Taxa_00104                      f__Microbacteriaceae
    ## Taxa_00105                      f__Microbacteriaceae
    ## Taxa_00106                      f__Microbacteriaceae
    ## Taxa_00107                      f__Microbacteriaceae
    ## Taxa_00108                      f__Microbacteriaceae
    ## Taxa_00109                      f__Microbacteriaceae
    ## Taxa_00110                      f__Microbacteriaceae
    ## Taxa_00111                      f__Microbacteriaceae
    ## Taxa_00112                      f__Microbacteriaceae
    ## Taxa_00113                      f__Microbacteriaceae
    ## Taxa_00114                      f__Microbacteriaceae
    ## Taxa_00115                      f__Microbacteriaceae
    ## Taxa_00116                         f__Micrococcaceae
    ## Taxa_00117                         f__Micrococcaceae
    ## Taxa_00118                         f__Micrococcaceae
    ## Taxa_00119                         f__Micrococcaceae
    ## Taxa_00120                         f__Micrococcaceae
    ## Taxa_00121                         f__Micrococcaceae
    ## Taxa_00122                     f__Micromonosporaceae
    ## Taxa_00123                     f__Micromonosporaceae
    ## Taxa_00124                     f__Micromonosporaceae
    ## Taxa_00125                     f__Micromonosporaceae
    ## Taxa_00126                     f__Micromonosporaceae
    ## Taxa_00127                     f__Micromonosporaceae
    ## Taxa_00128                     f__Micromonosporaceae
    ## Taxa_00129                     f__Micromonosporaceae
    ## Taxa_00130                     f__Micromonosporaceae
    ## Taxa_00131                     f__Micromonosporaceae
    ## Taxa_00132                     f__Micromonosporaceae
    ## Taxa_00133                     f__Micromonosporaceae
    ## Taxa_00134                     f__Micromonosporaceae
    ## Taxa_00135                     f__Micromonosporaceae
    ## Taxa_00136                     f__Micromonosporaceae
    ## Taxa_00137                     f__Micromonosporaceae
    ## Taxa_00138                     f__Micromonosporaceae
    ## Taxa_00139                     f__Micromonosporaceae
    ## Taxa_00140                     f__Micromonosporaceae
    ## Taxa_00141                     f__Micromonosporaceae
    ## Taxa_00142                     f__Micromonosporaceae
    ## Taxa_00143                     f__Micromonosporaceae
    ## Taxa_00144                     f__Micromonosporaceae
    ## Taxa_00145                     f__Micromonosporaceae
    ## Taxa_00146                     f__Micromonosporaceae
    ## Taxa_00147                       f__Mycobacteriaceae
    ## Taxa_00148                       f__Mycobacteriaceae
    ## Taxa_00149                        f__Nakamurellaceae
    ## Taxa_00150                           f__Nocardiaceae
    ## Taxa_00151                           f__Nocardiaceae
    ## Taxa_00152                           f__Nocardiaceae
    ## Taxa_00153                           f__Nocardiaceae
    ## Taxa_00154                           f__Nocardiaceae
    ## Taxa_00155                           f__Nocardiaceae
    ## Taxa_00156                           f__Nocardiaceae
    ## Taxa_00157                        f__Nocardioidaceae
    ## Taxa_00158                        f__Nocardioidaceae
    ## Taxa_00159                        f__Nocardioidaceae
    ## Taxa_00160                        f__Nocardioidaceae
    ## Taxa_00161                        f__Nocardioidaceae
    ## Taxa_00162                        f__Nocardioidaceae
    ## Taxa_00163                        f__Nocardioidaceae
    ## Taxa_00164                  f__Promicromonosporaceae
    ## Taxa_00165                  f__Promicromonosporaceae
    ## Taxa_00166                  f__Promicromonosporaceae
    ## Taxa_00167                   f__Propionibacteriaceae
    ## Taxa_00168                   f__Propionibacteriaceae
    ## Taxa_00169                   f__Propionibacteriaceae
    ## Taxa_00170                   f__Propionibacteriaceae
    ## Taxa_00171                   f__Propionibacteriaceae
    ## Taxa_00172                   f__Propionibacteriaceae
    ## Taxa_00173                     f__Pseudonocardiaceae
    ## Taxa_00174                     f__Pseudonocardiaceae
    ## Taxa_00175                     f__Pseudonocardiaceae
    ## Taxa_00176                     f__Pseudonocardiaceae
    ## Taxa_00177                     f__Pseudonocardiaceae
    ## Taxa_00178                     f__Pseudonocardiaceae
    ## Taxa_00179                     f__Pseudonocardiaceae
    ## Taxa_00180                     f__Pseudonocardiaceae
    ## Taxa_00181                     f__Pseudonocardiaceae
    ## Taxa_00182                     f__Pseudonocardiaceae
    ## Taxa_00183                     f__Pseudonocardiaceae
    ## Taxa_00184                     f__Pseudonocardiaceae
    ## Taxa_00185                     f__Pseudonocardiaceae
    ## Taxa_00186                     f__Pseudonocardiaceae
    ## Taxa_00187                     f__Pseudonocardiaceae
    ## Taxa_00188                     f__Pseudonocardiaceae
    ## Taxa_00189                     f__Pseudonocardiaceae
    ## Taxa_00190                     f__Pseudonocardiaceae
    ## Taxa_00191                             f__Ruaniaceae
    ## Taxa_00192                        f__Segniliparaceae
    ## Taxa_00193                        f__Sporichthyaceae
    ## Taxa_00194                      f__Streptomycetaceae
    ## Taxa_00195                      f__Streptomycetaceae
    ## Taxa_00196                      f__Streptomycetaceae
    ## Taxa_00197                      f__Streptomycetaceae
    ## Taxa_00198                   f__Streptosporangiaceae
    ## Taxa_00199                   f__Streptosporangiaceae
    ## Taxa_00200                   f__Streptosporangiaceae
    ## Taxa_00201                   f__Streptosporangiaceae
    ## Taxa_00202                   f__Streptosporangiaceae
    ## Taxa_00203                   f__Streptosporangiaceae
    ## Taxa_00204                   f__Streptosporangiaceae
    ## Taxa_00205     f__Streptosporangineae_incertae_sedis
    ## Taxa_00206                    f__Thermomonosporaceae
    ## Taxa_00207                    f__Thermomonosporaceae
    ## Taxa_00208                    f__Thermomonosporaceae
    ## Taxa_00209                    f__Thermomonosporaceae
    ## Taxa_00210                    f__Thermomonosporaceae
    ## Taxa_00211                    f__Thermomonosporaceae
    ## Taxa_00212                            f__Gaiellaceae
    ## Taxa_00213                                      <NA>
    ## Taxa_00214                      f__Conexibacteraceae
    ## Taxa_00215                      f__Patulibacteraceae
    ## Taxa_00216                   f__Solirubrobacteraceae
    ## Taxa_00217                     f__Thermoleophilaceae
    ## Taxa_00218                                      <NA>
    ## Taxa_00219                                      <NA>
    ## Taxa_00220             f__Aquificales_incertae_sedis
    ## Taxa_00221                                      <NA>
    ## Taxa_00222                    f__Armatimonadetes_gp2
    ## Taxa_00223                    f__Armatimonadetes_gp4
    ## Taxa_00224                    f__Armatimonadetes_gp5
    ## Taxa_00225                       f__Armatimonadaceae
    ## Taxa_00226                      f__Chthonomonadaceae
    ## Taxa_00227                      f__Fimbriimonadaceae
    ## Taxa_00228                                      <NA>
    ## Taxa_00229                                      <NA>
    ## Taxa_00230                      f__Marinilabiliaceae
    ## Taxa_00231                     f__Porphyromonadaceae
    ## Taxa_00232                     f__Porphyromonadaceae
    ## Taxa_00233                     f__Prolixibacteraceae
    ## Taxa_00234                                      <NA>
    ## Taxa_00235                           f__Chryseolinea
    ## Taxa_00236                      f__Cyclobacteriaceae
    ## Taxa_00237                      f__Cyclobacteriaceae
    ## Taxa_00238                      f__Cyclobacteriaceae
    ## Taxa_00239                          f__Cytophagaceae
    ## Taxa_00240                          f__Cytophagaceae
    ## Taxa_00241                          f__Cytophagaceae
    ## Taxa_00242                          f__Cytophagaceae
    ## Taxa_00243                          f__Cytophagaceae
    ## Taxa_00244                          f__Cytophagaceae
    ## Taxa_00245                          f__Cytophagaceae
    ## Taxa_00246                          f__Cytophagaceae
    ## Taxa_00247                          f__Cytophagaceae
    ## Taxa_00248                          f__Cytophagaceae
    ## Taxa_00249                          f__Cytophagaceae
    ## Taxa_00250                       f__Flammeovirgaceae
    ## Taxa_00251                       f__Flammeovirgaceae
    ## Taxa_00252                       f__Flammeovirgaceae
    ## Taxa_00253                       f__Flammeovirgaceae
    ## Taxa_00254                       f__Flammeovirgaceae
    ## Taxa_00255                       f__Flammeovirgaceae
    ## Taxa_00256                           f__Ohtaekwangia
    ## Taxa_00257                                      <NA>
    ## Taxa_00258                         f__Cryomorphaceae
    ## Taxa_00259                         f__Cryomorphaceae
    ## Taxa_00260                         f__Cryomorphaceae
    ## Taxa_00261                         f__Cryomorphaceae
    ## Taxa_00262                         f__Cryomorphaceae
    ## Taxa_00263                      f__Flavobacteriaceae
    ## Taxa_00264                      f__Flavobacteriaceae
    ## Taxa_00265                      f__Flavobacteriaceae
    ## Taxa_00266                      f__Flavobacteriaceae
    ## Taxa_00267                      f__Flavobacteriaceae
    ## Taxa_00268                      f__Flavobacteriaceae
    ## Taxa_00269                      f__Flavobacteriaceae
    ## Taxa_00270                                      <NA>
    ## Taxa_00271                       f__Chitinophagaceae
    ## Taxa_00272                       f__Chitinophagaceae
    ## Taxa_00273                       f__Chitinophagaceae
    ## Taxa_00274                       f__Chitinophagaceae
    ## Taxa_00275                       f__Chitinophagaceae
    ## Taxa_00276                       f__Chitinophagaceae
    ## Taxa_00277                       f__Chitinophagaceae
    ## Taxa_00278                       f__Chitinophagaceae
    ## Taxa_00279                       f__Chitinophagaceae
    ## Taxa_00280                       f__Chitinophagaceae
    ## Taxa_00281                       f__Chitinophagaceae
    ## Taxa_00282                       f__Chitinophagaceae
    ## Taxa_00283                       f__Chitinophagaceae
    ## Taxa_00284                       f__Chitinophagaceae
    ## Taxa_00285                       f__Chitinophagaceae
    ## Taxa_00286                       f__Chitinophagaceae
    ## Taxa_00287                       f__Chitinophagaceae
    ## Taxa_00288                       f__Chitinophagaceae
    ## Taxa_00289                       f__Chitinophagaceae
    ## Taxa_00290                       f__Chitinophagaceae
    ## Taxa_00291                       f__Chitinophagaceae
    ## Taxa_00292                       f__Chitinophagaceae
    ## Taxa_00293                         f__Saprospiraceae
    ## Taxa_00294                         f__Saprospiraceae
    ## Taxa_00295                         f__Saprospiraceae
    ## Taxa_00296                         f__Saprospiraceae
    ## Taxa_00297                    f__Sphingobacteriaceae
    ## Taxa_00298                    f__Sphingobacteriaceae
    ## Taxa_00299                    f__Sphingobacteriaceae
    ## Taxa_00300                    f__Sphingobacteriaceae
    ## Taxa_00301                    f__Sphingobacteriaceae
    ## Taxa_00302                    f__Sphingobacteriaceae
    ## Taxa_00303                    f__Sphingobacteriaceae
    ## Taxa_00304             f__BRC1_genera_incertae_sedis
    ## Taxa_00305            f__WPS-1_genera_incertae_sedis
    ## Taxa_00306            f__WPS-2_genera_incertae_sedis
    ## Taxa_00307              f__ZB3_genera_incertae_sedis
    ## Taxa_00308 f__Saccharibacteria_genera_incertae_sedis
    ## Taxa_00309                                      <NA>
    ## Taxa_00310                      f__Parachlamydiaceae
    ## Taxa_00311                      f__Parachlamydiaceae
    ## Taxa_00312                      f__Parachlamydiaceae
    ## Taxa_00313                           f__Simkaniaceae
    ## Taxa_00314                                      <NA>
    ## Taxa_00315                        f__Anaerolineaceae
    ## Taxa_00316                        f__Anaerolineaceae
    ## Taxa_00317                        f__Anaerolineaceae
    ## Taxa_00318                        f__Anaerolineaceae
    ## Taxa_00319                        f__Anaerolineaceae
    ## Taxa_00320                        f__Anaerolineaceae
    ## Taxa_00321                      f__Ardenticatenaceae
    ## Taxa_00322                         f__Caldilineaceae
    ## Taxa_00323                         f__Caldilineaceae
    ## Taxa_00324                         f__Caldilineaceae
    ## Taxa_00325                                      <NA>
    ## Taxa_00326                                      <NA>
    ## Taxa_00327                        f__Chloroflexaceae
    ## Taxa_00328                        f__Chloroflexaceae
    ## Taxa_00329                        f__Chloroflexaceae
    ## Taxa_00330                         f__Kallotenuaceae
    ## Taxa_00331                     f__Dehalococcoidaceae
    ## Taxa_00332                                      <NA>
    ## Taxa_00333                                      <NA>
    ## Taxa_00334                     f__Ktedonobacteraceae
    ## Taxa_00335                  f__Thermosporotrichaceae
    ## Taxa_00336                        f__Thermoflexaceae
    ## Taxa_00337                                      <NA>
    ## Taxa_00338                     f__Sphaerobacteraceae
    ## Taxa_00339                     f__Sphaerobacteraceae
    ## Taxa_00340                     f__Thermomicrobiaceae
    ## Taxa_00341                                      <NA>
    ## Taxa_00342                            f__Chloroplast
    ## Taxa_00343                            f__Chloroplast
    ## Taxa_00344                            f__Chloroplast
    ## Taxa_00345                            f__Chloroplast
    ## Taxa_00346                            f__Chloroplast
    ## Taxa_00347                            f__Chloroplast
    ## Taxa_00348                               f__Family I
    ## Taxa_00349                                      <NA>
    ## Taxa_00350                                      <NA>
    ## Taxa_00351                      f__Elusimicrobiaceae
    ## Taxa_00352               f__Candidatus Endomicrobium
    ## Taxa_00353                                      <NA>
    ## Taxa_00354                                      <NA>
    ## Taxa_00355                                      <NA>
    ## Taxa_00356                    f__Alicyclobacillaceae
    ## Taxa_00357                    f__Alicyclobacillaceae
    ## Taxa_00358                    f__Alicyclobacillaceae
    ## Taxa_00359                    f__Alicyclobacillaceae
    ## Taxa_00360                    f__Alicyclobacillaceae
    ## Taxa_00361                          f__Bacillaceae 1
    ## Taxa_00362                          f__Bacillaceae 1
    ## Taxa_00363                          f__Bacillaceae 1
    ## Taxa_00364                          f__Bacillaceae 1
    ## Taxa_00365                          f__Bacillaceae 1
    ## Taxa_00366                          f__Bacillaceae 1
    ## Taxa_00367                          f__Bacillaceae 1
    ## Taxa_00368                          f__Bacillaceae 2
    ## Taxa_00369                          f__Bacillaceae 2
    ## Taxa_00370                          f__Bacillaceae 2
    ## Taxa_00371                          f__Bacillaceae 2
    ## Taxa_00372                          f__Bacillaceae 2
    ## Taxa_00373                          f__Bacillaceae 2
    ## Taxa_00374                          f__Bacillaceae 2
    ## Taxa_00375                          f__Bacillaceae 2
    ## Taxa_00376                          f__Bacillaceae 2
    ## Taxa_00377                          f__Bacillaceae 2
    ## Taxa_00378                          f__Bacillaceae 2
    ## Taxa_00379                          f__Bacillaceae 2
    ## Taxa_00380                          f__Bacillaceae 2
    ## Taxa_00381            f__Bacillales_Incertae Sedis X
    ## Taxa_00382              f__Bacillales_incertae_sedis
    ## Taxa_00383              f__Bacillales_incertae_sedis
    ## Taxa_00384                     f__Paenibacillaceae 1
    ## Taxa_00385                     f__Paenibacillaceae 1
    ## Taxa_00386                     f__Paenibacillaceae 1
    ## Taxa_00387                     f__Paenibacillaceae 1
    ## Taxa_00388                     f__Paenibacillaceae 1
    ## Taxa_00389                     f__Paenibacillaceae 1
    ## Taxa_00390                     f__Paenibacillaceae 1
    ## Taxa_00391                     f__Paenibacillaceae 2
    ## Taxa_00392                     f__Paenibacillaceae 2
    ## Taxa_00393                     f__Paenibacillaceae 2
    ## Taxa_00394                     f__Paenibacillaceae 2
    ## Taxa_00395                         f__Planococcaceae
    ## Taxa_00396                         f__Planococcaceae
    ## Taxa_00397                         f__Planococcaceae
    ## Taxa_00398                         f__Planococcaceae
    ## Taxa_00399                         f__Planococcaceae
    ## Taxa_00400                         f__Planococcaceae
    ## Taxa_00401                         f__Planococcaceae
    ## Taxa_00402                         f__Planococcaceae
    ## Taxa_00403                         f__Planococcaceae
    ## Taxa_00404                         f__Planococcaceae
    ## Taxa_00405                         f__Planococcaceae
    ## Taxa_00406                         f__Planococcaceae
    ## Taxa_00407                         f__Planococcaceae
    ## Taxa_00408                      f__Staphylococcaceae
    ## Taxa_00409               f__Thermoactinomycetaceae 1
    ## Taxa_00410               f__Thermoactinomycetaceae 1
    ## Taxa_00411               f__Thermoactinomycetaceae 1
    ## Taxa_00412               f__Thermoactinomycetaceae 1
    ## Taxa_00413               f__Thermoactinomycetaceae 1
    ## Taxa_00414               f__Thermoactinomycetaceae 1
    ## Taxa_00415               f__Thermoactinomycetaceae 2
    ## Taxa_00416               f__Thermoactinomycetaceae 2
    ## Taxa_00417               f__Thermoactinomycetaceae 2
    ## Taxa_00418                                      <NA>
    ## Taxa_00419                      f__Carnobacteriaceae
    ## Taxa_00420                        f__Enterococcaceae
    ## Taxa_00421                       f__Lactobacillaceae
    ## Taxa_00422                       f__Streptococcaceae
    ## Taxa_00423                       f__Streptococcaceae
    ## Taxa_00424                                      <NA>
    ## Taxa_00425                                      <NA>
    ## Taxa_00426                       f__Catabacteriaceae
    ## Taxa_00427                       f__Clostridiaceae 1
    ## Taxa_00428                       f__Clostridiaceae 1
    ## Taxa_00429                       f__Clostridiaceae 1
    ## Taxa_00430                       f__Clostridiaceae 1
    ## Taxa_00431                       f__Clostridiaceae 1
    ## Taxa_00432                       f__Clostridiaceae 1
    ## Taxa_00433                       f__Clostridiaceae 1
    ## Taxa_00434                       f__Clostridiaceae 4
    ## Taxa_00435       f__Clostridiales_Incertae Sedis III
    ## Taxa_00436       f__Clostridiales_Incertae Sedis III
    ## Taxa_00437        f__Clostridiales_Incertae Sedis XI
    ## Taxa_00438      f__Clostridiales_Incertae Sedis XIII
    ## Taxa_00439     f__Clostridiales_Incertae Sedis XVIII
    ## Taxa_00440                     f__Gracilibacteraceae
    ## Taxa_00441                     f__Gracilibacteraceae
    ## Taxa_00442                      f__Heliobacteriaceae
    ## Taxa_00443                        f__Lachnospiraceae
    ## Taxa_00444                        f__Lachnospiraceae
    ## Taxa_00445                        f__Lachnospiraceae
    ## Taxa_00446                        f__Lachnospiraceae
    ## Taxa_00447                       f__Peptococcaceae 1
    ## Taxa_00448                       f__Peptococcaceae 1
    ## Taxa_00449                        f__Ruminococcaceae
    ## Taxa_00450                        f__Ruminococcaceae
    ## Taxa_00451                        f__Ruminococcaceae
    ## Taxa_00452                        f__Ruminococcaceae
    ## Taxa_00453                        f__Ruminococcaceae
    ## Taxa_00454                        f__Ruminococcaceae
    ## Taxa_00455                        f__Ruminococcaceae
    ## Taxa_00456                        f__Ruminococcaceae
    ## Taxa_00457                        f__Ruminococcaceae
    ## Taxa_00458                        f__Ruminococcaceae
    ## Taxa_00459                    f__Syntrophomonadaceae
    ## Taxa_00460                   f__Thermodesulfobiaceae
    ## Taxa_00461                    f__Erysipelotrichaceae
    ## Taxa_00462                                      <NA>
    ## Taxa_00463                     f__Acidaminococcaceae
    ## Taxa_00464                        f__Veillonellaceae
    ## Taxa_00465                        f__Veillonellaceae
    ## Taxa_00466                        f__Veillonellaceae
    ## Taxa_00467                      f__Gemmatimonadaceae
    ## Taxa_00468               f__Candidatus Hydrogenedens
    ## Taxa_00469  f__Latescibacteria_genera_incertae_sedis
    ## Taxa_00470                                      <NA>
    ## Taxa_00471                         f__Victivallaceae
    ## Taxa_00472   f__Microgenomates_genera_incertae_sedis
    ## Taxa_00473                         f__Nitrospiraceae
    ## Taxa_00474     f__Omnitrophica_genera_incertae_sedis
    ## Taxa_00475    f__Parcubacteria_genera_incertae_sedis
    ## Taxa_00476                                      <NA>
    ## Taxa_00477                                      <NA>
    ## Taxa_00478                       f__Phycisphaeraceae
    ## Taxa_00479                      f__Tepidisphaeraceae
    ## Taxa_00480                                      <NA>
    ## Taxa_00481                f__Candidatus Brocadiaceae
    ## Taxa_00482                f__Candidatus Brocadiaceae
    ## Taxa_00483                      f__Planctomycetaceae
    ## Taxa_00484                      f__Planctomycetaceae
    ## Taxa_00485                      f__Planctomycetaceae
    ## Taxa_00486                      f__Planctomycetaceae
    ## Taxa_00487                      f__Planctomycetaceae
    ## Taxa_00488                      f__Planctomycetaceae
    ## Taxa_00489                      f__Planctomycetaceae
    ## Taxa_00490                      f__Planctomycetaceae
    ## Taxa_00491                      f__Planctomycetaceae
    ## Taxa_00492                      f__Planctomycetaceae
    ## Taxa_00493                      f__Planctomycetaceae
    ## Taxa_00494                      f__Planctomycetaceae
    ## Taxa_00495                      f__Planctomycetaceae
    ## Taxa_00496                      f__Planctomycetaceae
    ## Taxa_00497                      f__Planctomycetaceae
    ## Taxa_00498                                      <NA>
    ## Taxa_00499                                      <NA>
    ## Taxa_00500                                      <NA>
    ## Taxa_00501                             f__Breoghania
    ## Taxa_00502                           f__Geminicoccus
    ## Taxa_00503                         f__Rhizomicrobium
    ## Taxa_00504                                      <NA>
    ## Taxa_00505                       f__Caulobacteraceae
    ## Taxa_00506                       f__Caulobacteraceae
    ## Taxa_00507                       f__Caulobacteraceae
    ## Taxa_00508                       f__Caulobacteraceae
    ## Taxa_00509                       f__Caulobacteraceae
    ## Taxa_00510                        f__Hyphomonadaceae
    ## Taxa_00511                        f__Hyphomonadaceae
    ## Taxa_00512                        f__Hyphomonadaceae
    ## Taxa_00513                            f__Eilatimonas
    ## Taxa_00514                       f__Parvularculaceae
    ## Taxa_00515                                      <NA>
    ## Taxa_00516                      f__Aurantimonadaceae
    ## Taxa_00517                       f__Beijerinckiaceae
    ## Taxa_00518                       f__Beijerinckiaceae
    ## Taxa_00519                       f__Beijerinckiaceae
    ## Taxa_00520                       f__Beijerinckiaceae
    ## Taxa_00521                       f__Beijerinckiaceae
    ## Taxa_00522                       f__Beijerinckiaceae
    ## Taxa_00523                       f__Beijerinckiaceae
    ## Taxa_00524                       f__Beijerinckiaceae
    ## Taxa_00525                       f__Beijerinckiaceae
    ## Taxa_00526                       f__Beijerinckiaceae
    ## Taxa_00527                      f__Bradyrhizobiaceae
    ## Taxa_00528                      f__Bradyrhizobiaceae
    ## Taxa_00529                      f__Bradyrhizobiaceae
    ## Taxa_00530                      f__Bradyrhizobiaceae
    ## Taxa_00531                      f__Bradyrhizobiaceae
    ## Taxa_00532                      f__Bradyrhizobiaceae
    ## Taxa_00533                      f__Bradyrhizobiaceae
    ## Taxa_00534                      f__Bradyrhizobiaceae
    ## Taxa_00535                      f__Bradyrhizobiaceae
    ## Taxa_00536                      f__Bradyrhizobiaceae
    ## Taxa_00537                           f__Brucellaceae
    ## Taxa_00538                           f__Brucellaceae
    ## Taxa_00539                           f__Brucellaceae
    ## Taxa_00540                           f__Brucellaceae
    ## Taxa_00541                      f__Hyphomicrobiaceae
    ## Taxa_00542                      f__Hyphomicrobiaceae
    ## Taxa_00543                      f__Hyphomicrobiaceae
    ## Taxa_00544                      f__Hyphomicrobiaceae
    ## Taxa_00545                      f__Hyphomicrobiaceae
    ## Taxa_00546                      f__Hyphomicrobiaceae
    ## Taxa_00547                      f__Hyphomicrobiaceae
    ## Taxa_00548                      f__Hyphomicrobiaceae
    ## Taxa_00549                      f__Hyphomicrobiaceae
    ## Taxa_00550                      f__Hyphomicrobiaceae
    ## Taxa_00551                      f__Hyphomicrobiaceae
    ## Taxa_00552                      f__Hyphomicrobiaceae
    ## Taxa_00553                    f__Methylobacteriaceae
    ## Taxa_00554                    f__Methylobacteriaceae
    ## Taxa_00555                    f__Methylobacteriaceae
    ## Taxa_00556                    f__Methylobacteriaceae
    ## Taxa_00557                    f__Methylobacteriaceae
    ## Taxa_00558                       f__Methylocystaceae
    ## Taxa_00559                       f__Methylocystaceae
    ## Taxa_00560                       f__Methylocystaceae
    ## Taxa_00561                       f__Methylocystaceae
    ## Taxa_00562                       f__Methylocystaceae
    ## Taxa_00563                       f__Methylocystaceae
    ## Taxa_00564                     f__Phyllobacteriaceae
    ## Taxa_00565                     f__Phyllobacteriaceae
    ## Taxa_00566                     f__Phyllobacteriaceae
    ## Taxa_00567                     f__Phyllobacteriaceae
    ## Taxa_00568                     f__Phyllobacteriaceae
    ## Taxa_00569                     f__Phyllobacteriaceae
    ## Taxa_00570                     f__Phyllobacteriaceae
    ## Taxa_00571                     f__Phyllobacteriaceae
    ## Taxa_00572                     f__Phyllobacteriaceae
    ## Taxa_00573                           f__Rhizobiaceae
    ## Taxa_00574                           f__Rhizobiaceae
    ## Taxa_00575                           f__Rhizobiaceae
    ## Taxa_00576                           f__Rhizobiaceae
    ## Taxa_00577                           f__Rhizobiaceae
    ## Taxa_00578                           f__Rhizobiaceae
    ## Taxa_00579                           f__Rhizobiaceae
    ## Taxa_00580             f__Rhizobiales_incertae_sedis
    ## Taxa_00581             f__Rhizobiales_incertae_sedis
    ## Taxa_00582             f__Rhizobiales_incertae_sedis
    ## Taxa_00583             f__Rhizobiales_incertae_sedis
    ## Taxa_00584             f__Rhizobiales_incertae_sedis
    ## Taxa_00585             f__Rhizobiales_incertae_sedis
    ## Taxa_00586                           f__Rhodobiaceae
    ## Taxa_00587                           f__Rhodobiaceae
    ## Taxa_00588                           f__Rhodobiaceae
    ## Taxa_00589                           f__Rhodobiaceae
    ## Taxa_00590                           f__Rhodobiaceae
    ## Taxa_00591                           f__Rhodobiaceae
    ## Taxa_00592                           f__Rhodobiaceae
    ## Taxa_00593                           f__Rhodobiaceae
    ## Taxa_00594                           f__Rhodobiaceae
    ## Taxa_00595                          f__Roseiarcaceae
    ## Taxa_00596                      f__Xanthobacteraceae
    ## Taxa_00597                      f__Xanthobacteraceae
    ## Taxa_00598                      f__Xanthobacteraceae
    ## Taxa_00599                       f__Rhodobacteraceae
    ## Taxa_00600                       f__Rhodobacteraceae
    ## Taxa_00601                       f__Rhodobacteraceae
    ## Taxa_00602                       f__Rhodobacteraceae
    ## Taxa_00603                       f__Rhodobacteraceae
    ## Taxa_00604                       f__Rhodobacteraceae
    ## Taxa_00605                       f__Rhodobacteraceae
    ## Taxa_00606                       f__Rhodobacteraceae
    ## Taxa_00607                                      <NA>
    ## Taxa_00608                       f__Acetobacteraceae
    ## Taxa_00609                       f__Acetobacteraceae
    ## Taxa_00610                       f__Acetobacteraceae
    ## Taxa_00611                       f__Acetobacteraceae
    ## Taxa_00612                       f__Acetobacteraceae
    ## Taxa_00613                       f__Acetobacteraceae
    ## Taxa_00614                       f__Acetobacteraceae
    ## Taxa_00615                       f__Acetobacteraceae
    ## Taxa_00616                       f__Acetobacteraceae
    ## Taxa_00617                       f__Acetobacteraceae
    ## Taxa_00618                       f__Acetobacteraceae
    ## Taxa_00619                       f__Acetobacteraceae
    ## Taxa_00620                       f__Acetobacteraceae
    ## Taxa_00621                       f__Acetobacteraceae
    ## Taxa_00622                       f__Acetobacteraceae
    ## Taxa_00623                       f__Acetobacteraceae
    ## Taxa_00624                       f__Acetobacteraceae
    ## Taxa_00625                       f__Acetobacteraceae
    ## Taxa_00626                       f__Acetobacteraceae
    ## Taxa_00627                       f__Acetobacteraceae
    ## Taxa_00628                       f__Acetobacteraceae
    ## Taxa_00629                       f__Acetobacteraceae
    ## Taxa_00630                       f__Acetobacteraceae
    ## Taxa_00631                       f__Acetobacteraceae
    ## Taxa_00632                               f__Elioraea
    ## Taxa_00633                             f__Reyranella
    ## Taxa_00634                      f__Rhodospirillaceae
    ## Taxa_00635                      f__Rhodospirillaceae
    ## Taxa_00636                      f__Rhodospirillaceae
    ## Taxa_00637                      f__Rhodospirillaceae
    ## Taxa_00638                      f__Rhodospirillaceae
    ## Taxa_00639                      f__Rhodospirillaceae
    ## Taxa_00640                      f__Rhodospirillaceae
    ## Taxa_00641                      f__Rhodospirillaceae
    ## Taxa_00642                      f__Rhodospirillaceae
    ## Taxa_00643                      f__Rhodospirillaceae
    ## Taxa_00644                      f__Rhodospirillaceae
    ## Taxa_00645                      f__Rhodospirillaceae
    ## Taxa_00646                      f__Rhodospirillaceae
    ## Taxa_00647                      f__Rhodospirillaceae
    ## Taxa_00648                      f__Rhodospirillaceae
    ## Taxa_00649                      f__Rhodospirillaceae
    ## Taxa_00650                      f__Rhodospirillaceae
    ## Taxa_00651                      f__Rhodospirillaceae
    ## Taxa_00652                      f__Rhodospirillaceae
    ## Taxa_00653                                      <NA>
    ## Taxa_00654                        f__Anaplasmataceae
    ## Taxa_00655                         f__Rickettsiaceae
    ## Taxa_00656                         f__Rickettsiaceae
    ## Taxa_00657                         f__Rickettsiaceae
    ## Taxa_00658                        f__Sneathiellaceae
    ## Taxa_00659                        f__Sneathiellaceae
    ## Taxa_00660                                      <NA>
    ## Taxa_00661                     f__Erythrobacteraceae
    ## Taxa_00662                     f__Erythrobacteraceae
    ## Taxa_00663                     f__Erythrobacteraceae
    ## Taxa_00664                     f__Erythrobacteraceae
    ## Taxa_00665                      f__Sphingomonadaceae
    ## Taxa_00666                      f__Sphingomonadaceae
    ## Taxa_00667                      f__Sphingomonadaceae
    ## Taxa_00668                      f__Sphingomonadaceae
    ## Taxa_00669                      f__Sphingomonadaceae
    ## Taxa_00670                      f__Sphingomonadaceae
    ## Taxa_00671                      f__Sphingomonadaceae
    ## Taxa_00672                      f__Sphingomonadaceae
    ## Taxa_00673                      f__Sphingomonadaceae
    ## Taxa_00674                      f__Sphingomonadaceae
    ## Taxa_00675                      f__Sphingomonadaceae
    ## Taxa_00676                      f__Sphingomonadaceae
    ## Taxa_00677                      f__Sphingomonadaceae
    ## Taxa_00678                      f__Sphingomonadaceae
    ## Taxa_00679                      f__Sphingomonadaceae
    ## Taxa_00680                      f__Sphingomonadaceae
    ## Taxa_00681                      f__Sphingomonadaceae
    ## Taxa_00682                                      <NA>
    ## Taxa_00683                                      <NA>
    ## Taxa_00684                         f__Alcaligenaceae
    ## Taxa_00685                         f__Alcaligenaceae
    ## Taxa_00686                         f__Alcaligenaceae
    ## Taxa_00687                         f__Alcaligenaceae
    ## Taxa_00688                         f__Alcaligenaceae
    ## Taxa_00689                         f__Alcaligenaceae
    ## Taxa_00690                         f__Alcaligenaceae
    ## Taxa_00691                         f__Alcaligenaceae
    ## Taxa_00692                         f__Alcaligenaceae
    ## Taxa_00693                         f__Alcaligenaceae
    ## Taxa_00694                       f__Burkholderiaceae
    ## Taxa_00695                       f__Burkholderiaceae
    ## Taxa_00696                       f__Burkholderiaceae
    ## Taxa_00697                       f__Burkholderiaceae
    ## Taxa_00698                       f__Burkholderiaceae
    ## Taxa_00699                       f__Burkholderiaceae
    ## Taxa_00700         f__Burkholderiales_incertae_sedis
    ## Taxa_00701         f__Burkholderiales_incertae_sedis
    ## Taxa_00702         f__Burkholderiales_incertae_sedis
    ## Taxa_00703         f__Burkholderiales_incertae_sedis
    ## Taxa_00704         f__Burkholderiales_incertae_sedis
    ## Taxa_00705         f__Burkholderiales_incertae_sedis
    ## Taxa_00706         f__Burkholderiales_incertae_sedis
    ## Taxa_00707         f__Burkholderiales_incertae_sedis
    ## Taxa_00708         f__Burkholderiales_incertae_sedis
    ## Taxa_00709         f__Burkholderiales_incertae_sedis
    ## Taxa_00710         f__Burkholderiales_incertae_sedis
    ## Taxa_00711         f__Burkholderiales_incertae_sedis
    ## Taxa_00712                         f__Comamonadaceae
    ## Taxa_00713                         f__Comamonadaceae
    ## Taxa_00714                         f__Comamonadaceae
    ## Taxa_00715                         f__Comamonadaceae
    ## Taxa_00716                         f__Comamonadaceae
    ## Taxa_00717                         f__Comamonadaceae
    ## Taxa_00718                         f__Comamonadaceae
    ## Taxa_00719                         f__Comamonadaceae
    ## Taxa_00720                         f__Comamonadaceae
    ## Taxa_00721                         f__Comamonadaceae
    ## Taxa_00722                         f__Comamonadaceae
    ## Taxa_00723                         f__Comamonadaceae
    ## Taxa_00724                         f__Comamonadaceae
    ## Taxa_00725                         f__Comamonadaceae
    ## Taxa_00726                         f__Comamonadaceae
    ## Taxa_00727                         f__Comamonadaceae
    ## Taxa_00728                         f__Comamonadaceae
    ## Taxa_00729                         f__Comamonadaceae
    ## Taxa_00730                         f__Comamonadaceae
    ## Taxa_00731                         f__Comamonadaceae
    ## Taxa_00732                         f__Comamonadaceae
    ## Taxa_00733                         f__Comamonadaceae
    ## Taxa_00734                         f__Comamonadaceae
    ## Taxa_00735                       f__Oxalobacteraceae
    ## Taxa_00736                       f__Oxalobacteraceae
    ## Taxa_00737                       f__Oxalobacteraceae
    ## Taxa_00738                       f__Oxalobacteraceae
    ## Taxa_00739                       f__Oxalobacteraceae
    ## Taxa_00740                       f__Oxalobacteraceae
    ## Taxa_00741                       f__Oxalobacteraceae
    ## Taxa_00742                       f__Oxalobacteraceae
    ## Taxa_00743                       f__Oxalobacteraceae
    ## Taxa_00744                       f__Oxalobacteraceae
    ## Taxa_00745                       f__Oxalobacteraceae
    ## Taxa_00746                       f__Oxalobacteraceae
    ## Taxa_00747                       f__Oxalobacteraceae
    ## Taxa_00748                       f__Oxalobacteraceae
    ## Taxa_00749                         f__Sutterellaceae
    ## Taxa_00750                      f__Ferritrophicaceae
    ## Taxa_00751                            f__Ferrovaceae
    ## Taxa_00752                        f__Gallionellaceae
    ## Taxa_00753                        f__Gallionellaceae
    ## Taxa_00754                        f__Gallionellaceae
    ## Taxa_00755                       f__Methylophilaceae
    ## Taxa_00756                       f__Methylophilaceae
    ## Taxa_00757                          f__Neisseriaceae
    ## Taxa_00758                          f__Neisseriaceae
    ## Taxa_00759                          f__Neisseriaceae
    ## Taxa_00760                          f__Neisseriaceae
    ## Taxa_00761                          f__Neisseriaceae
    ## Taxa_00762                          f__Neisseriaceae
    ## Taxa_00763                          f__Neisseriaceae
    ## Taxa_00764                      f__Nitrosomonadaceae
    ## Taxa_00765                      f__Nitrosomonadaceae
    ## Taxa_00766                      f__Procabacteriaceae
    ## Taxa_00767                         f__Rhodocyclaceae
    ## Taxa_00768                         f__Rhodocyclaceae
    ## Taxa_00769                         f__Rhodocyclaceae
    ## Taxa_00770                         f__Rhodocyclaceae
    ## Taxa_00771                         f__Rhodocyclaceae
    ## Taxa_00772                         f__Rhodocyclaceae
    ## Taxa_00773                         f__Rhodocyclaceae
    ## Taxa_00774                         f__Rhodocyclaceae
    ## Taxa_00775                         f__Rhodocyclaceae
    ## Taxa_00776                         f__Rhodocyclaceae
    ## Taxa_00777                       f__Sulfuricellaceae
    ## Taxa_00778                       f__Sulfuricellaceae
    ## Taxa_00779                                      <NA>
    ## Taxa_00780                                      <NA>
    ## Taxa_00781                     f__Bacteriovoracaceae
    ## Taxa_00782                     f__Bacteriovoracaceae
    ## Taxa_00783                     f__Bacteriovoracaceae
    ## Taxa_00784                     f__Bdellovibrionaceae
    ## Taxa_00785                     f__Bdellovibrionaceae
    ## Taxa_00786                     f__Bdellovibrionaceae
    ## Taxa_00787               f__Pseudobacteriovoracaceae
    ## Taxa_00788                                      <NA>
    ## Taxa_00789                            f__Deferrisoma
    ## Taxa_00790                       f__Dissulfuribacter
    ## Taxa_00791                                      <NA>
    ## Taxa_00792                     f__Desulfobacteraceae
    ## Taxa_00793                     f__Desulfobacteraceae
    ## Taxa_00794                    f__Desulfovibrionaceae
    ## Taxa_00795                    f__Desulfovibrionaceae
    ## Taxa_00796                                      <NA>
    ## Taxa_00797                         f__Geobacteraceae
    ## Taxa_00798                         f__Geobacteraceae
    ## Taxa_00799                         f__Geobacteraceae
    ## Taxa_00800                                      <NA>
    ## Taxa_00801                       f__Cystobacteraceae
    ## Taxa_00802                       f__Cystobacteraceae
    ## Taxa_00803                       f__Cystobacteraceae
    ## Taxa_00804                       f__Cystobacteraceae
    ## Taxa_00805                       f__Cystobacteraceae
    ## Taxa_00806                       f__Cystobacteraceae
    ## Taxa_00807                       f__Cystobacteraceae
    ## Taxa_00808                          f__Haliangiaceae
    ## Taxa_00809                           f__Kofleriaceae
    ## Taxa_00810                       f__Labilitrichaceae
    ## Taxa_00811                          f__Myxococcaceae
    ## Taxa_00812                          f__Myxococcaceae
    ## Taxa_00813                          f__Myxococcaceae
    ## Taxa_00814                          f__Myxococcaceae
    ## Taxa_00815                         f__Nannocystaceae
    ## Taxa_00816                         f__Nannocystaceae
    ## Taxa_00817                         f__Nannocystaceae
    ## Taxa_00818                     f__Phaselicystidaceae
    ## Taxa_00819                          f__Polyangiaceae
    ## Taxa_00820                          f__Polyangiaceae
    ## Taxa_00821                          f__Polyangiaceae
    ## Taxa_00822                          f__Polyangiaceae
    ## Taxa_00823                          f__Polyangiaceae
    ## Taxa_00824                          f__Polyangiaceae
    ## Taxa_00825                        f__Sandaracinaceae
    ## Taxa_00826                     f__Vulgatibacteraceae
    ## Taxa_00827                                      <NA>
    ## Taxa_00828                          f__Syntrophaceae
    ## Taxa_00829                   f__Syntrophobacteraceae
    ## Taxa_00830                                      <NA>
    ## Taxa_00831                           f__Nautiliaceae
    ## Taxa_00832                                      <NA>
    ## Taxa_00833                                      <NA>
    ## Taxa_00834                          f__Chromatiaceae
    ## Taxa_00835                          f__Chromatiaceae
    ## Taxa_00836                          f__Chromatiaceae
    ## Taxa_00837                          f__Chromatiaceae
    ## Taxa_00838                          f__Chromatiaceae
    ## Taxa_00839                          f__Chromatiaceae
    ## Taxa_00840                 f__Ectothiorhodospiraceae
    ## Taxa_00841                 f__Ectothiorhodospiraceae
    ## Taxa_00842                 f__Ectothiorhodospiraceae
    ## Taxa_00843                 f__Ectothiorhodospiraceae
    ## Taxa_00844                 f__Ectothiorhodospiraceae
    ## Taxa_00845                 f__Ectothiorhodospiraceae
    ## Taxa_00846                 f__Ectothiorhodospiraceae
    ## Taxa_00847                 f__Ectothiorhodospiraceae
    ## Taxa_00848                     f__Enterobacteriaceae
    ## Taxa_00849                     f__Enterobacteriaceae
    ## Taxa_00850                     f__Enterobacteriaceae
    ## Taxa_00851                     f__Enterobacteriaceae
    ## Taxa_00852                     f__Enterobacteriaceae
    ## Taxa_00853                                      <NA>
    ## Taxa_00854                  f__Candidatus Carsonella
    ## Taxa_00855                       f__Methylohalomonas
    ## Taxa_00856                          f__Sedimenticola
    ## Taxa_00857                                      <NA>
    ## Taxa_00858                           f__Coxiellaceae
    ## Taxa_00859                           f__Coxiellaceae
    ## Taxa_00860                           f__Coxiellaceae
    ## Taxa_00861                           f__Coxiellaceae
    ## Taxa_00862                         f__Legionellaceae
    ## Taxa_00863                       f__Methylococcaceae
    ## Taxa_00864                       f__Methylococcaceae
    ## Taxa_00865                       f__Methylococcaceae
    ## Taxa_00866                       f__Methylococcaceae
    ## Taxa_00867                       f__Methylococcaceae
    ## Taxa_00868                                      <NA>
    ## Taxa_00869                            f__Hahellaceae
    ## Taxa_00870                     f__Oceanospirillaceae
    ## Taxa_00871       f__Oceanospirillales_incertae_sedis
    ## Taxa_00872                        f__Pasteurellaceae
    ## Taxa_00873                                      <NA>
    ## Taxa_00874                          f__Moraxellaceae
    ## Taxa_00875                       f__Pseudomonadaceae
    ## Taxa_00876                       f__Pseudomonadaceae
    ## Taxa_00877                       f__Pseudomonadaceae
    ## Taxa_00878                       f__Pseudomonadaceae
    ## Taxa_00879                       f__Pseudomonadaceae
    ## Taxa_00880                       f__Pseudomonadaceae
    ## Taxa_00881                                      <NA>
    ## Taxa_00882                         f__Thiotrichaceae
    ## Taxa_00883                                      <NA>
    ## Taxa_00884                          f__Algiphilaceae
    ## Taxa_00885                        f__Sinobacteraceae
    ## Taxa_00886                        f__Sinobacteraceae
    ## Taxa_00887                        f__Sinobacteraceae
    ## Taxa_00888                        f__Sinobacteraceae
    ## Taxa_00889                        f__Sinobacteraceae
    ## Taxa_00890                        f__Sinobacteraceae
    ## Taxa_00891                        f__Sinobacteraceae
    ## Taxa_00892                        f__Sinobacteraceae
    ## Taxa_00893                       f__Xanthomonadaceae
    ## Taxa_00894                       f__Xanthomonadaceae
    ## Taxa_00895                       f__Xanthomonadaceae
    ## Taxa_00896                       f__Xanthomonadaceae
    ## Taxa_00897                       f__Xanthomonadaceae
    ## Taxa_00898                       f__Xanthomonadaceae
    ## Taxa_00899                       f__Xanthomonadaceae
    ## Taxa_00900                       f__Xanthomonadaceae
    ## Taxa_00901                       f__Xanthomonadaceae
    ## Taxa_00902                       f__Xanthomonadaceae
    ## Taxa_00903                       f__Xanthomonadaceae
    ## Taxa_00904                       f__Xanthomonadaceae
    ## Taxa_00905                       f__Xanthomonadaceae
    ## Taxa_00906                       f__Xanthomonadaceae
    ## Taxa_00907                       f__Xanthomonadaceae
    ## Taxa_00908                       f__Xanthomonadaceae
    ## Taxa_00909                       f__Xanthomonadaceae
    ## Taxa_00910                       f__Xanthomonadaceae
    ## Taxa_00911                       f__Xanthomonadaceae
    ## Taxa_00912                       f__Xanthomonadaceae
    ## Taxa_00913                       f__Xanthomonadaceae
    ## Taxa_00914                         f__Oligoflexaceae
    ## Taxa_00915                                      <NA>
    ## Taxa_00916                        f__Brevinemataceae
    ## Taxa_00917                         f__Leptospiraceae
    ## Taxa_00918                         f__Leptospiraceae
    ## Taxa_00919                         f__Leptospiraceae
    ## Taxa_00920                        f__Spirochaetaceae
    ## Taxa_00921                        f__Spirochaetaceae
    ## Taxa_00922                        f__Spirochaetaceae
    ## Taxa_00923              f__SR1_genera_incertae_sedis
    ## Taxa_00924                     f__Entomoplasmataceae
    ## Taxa_00925                       f__Haloplasmataceae
    ## Taxa_00926              f__Thermodesulfobacteriaceae
    ## Taxa_00927                                      <NA>
    ## Taxa_00928                                      <NA>
    ## Taxa_00929                            f__Opitutaceae
    ## Taxa_00930                            f__Opitutaceae
    ## Taxa_00931                            f__Opitutaceae
    ## Taxa_00932                       f__Puniceicoccaceae
    ## Taxa_00933                       f__Puniceicoccaceae
    ## Taxa_00934                       f__Puniceicoccaceae
    ## Taxa_00935                                      <NA>
    ## Taxa_00936   f__Spartobacteria_genera_incertae_sedis
    ## Taxa_00937                         f__Terrimicrobium
    ## Taxa_00938                                      <NA>
    ## Taxa_00939                            f__Limisphaera
    ## Taxa_00940     f__Subdivision3_genera_incertae_sedis
    ## Taxa_00941                    f__Verrucomicrobiaceae
    ## Taxa_00942                    f__Verrucomicrobiaceae
    ## Taxa_00943                    f__Verrucomicrobiaceae
    ## Taxa_00944                    f__Verrucomicrobiaceae
    ## Taxa_00945                    f__Verrucomicrobiaceae
    ## Taxa_00946                    f__Verrucomicrobiaceae
    ## Taxa_00947                    f__Verrucomicrobiaceae
    ## Taxa_00948                    f__Verrucomicrobiaceae
    ## Taxa_00949                                      <NA>
    ##                                                Genus
    ## Taxa_00000                                      <NA>
    ## Taxa_00001                              g__Pyrolobus
    ## Taxa_00002                                      <NA>
    ## Taxa_00003                                      <NA>
    ## Taxa_00004                                      <NA>
    ## Taxa_00005                                   g__Gp10
    ## Taxa_00006                                   g__Gp11
    ## Taxa_00007                                   g__Gp12
    ## Taxa_00008                                   g__Gp13
    ## Taxa_00009                                   g__Gp15
    ## Taxa_00010                                   g__Gp16
    ## Taxa_00011                                   g__Gp17
    ## Taxa_00012                                   g__Gp18
    ## Taxa_00013                                   g__Gp19
    ## Taxa_00014                             g__Acidicapsa
    ## Taxa_00015                              g__Acidipila
    ## Taxa_00016                         g__Acidobacterium
    ## Taxa_00017                              g__Bryocella
    ## Taxa_00018                  g__Candidatus Koribacter
    ## Taxa_00019                           g__Edaphobacter
    ## Taxa_00020                                    g__Gp1
    ## Taxa_00021                           g__Granulicella
    ## Taxa_00022                          g__Telmatobacter
    ## Taxa_00023                            g__Terriglobus
    ## Taxa_00024                                   g__Gp20
    ## Taxa_00025                                   g__Gp22
    ## Taxa_00026                                   g__Gp25
    ## Taxa_00027                                    g__Gp2
    ## Taxa_00028                                      <NA>
    ## Taxa_00029                             g__Bryobacter
    ## Taxa_00030                  g__Candidatus Solibacter
    ## Taxa_00031                                    g__Gp3
    ## Taxa_00032                          g__Paludibaculum
    ## Taxa_00033                                      <NA>
    ## Taxa_00034                            g__Aridibacter
    ## Taxa_00035                          g__Blastocatella
    ## Taxa_00036                                    g__Gp4
    ## Taxa_00037                            g__Pyrinomonas
    ## Taxa_00038                                    g__Gp5
    ## Taxa_00039                                    g__Gp6
    ## Taxa_00040                                    g__Gp7
    ## Taxa_00041                                      <NA>
    ## Taxa_00042                               g__Geothrix
    ## Taxa_00043                              g__Holophaga
    ## Taxa_00044                                      <NA>
    ## Taxa_00045                                      <NA>
    ## Taxa_00046                                      <NA>
    ## Taxa_00047                                      <NA>
    ## Taxa_00048                         g__Ferrimicrobium
    ## Taxa_00049                             g__Ferrithrix
    ## Taxa_00050                          g__Ilumatobacter
    ## Taxa_00051                        g__Aciditerrimonas
    ## Taxa_00052                                      <NA>
    ## Taxa_00053                           g__Aquihabitans
    ## Taxa_00054                                  g__Iamia
    ## Taxa_00055                                      <NA>
    ## Taxa_00056                           g__Acidothermus
    ## Taxa_00057                            g__Actinospica
    ## Taxa_00058                                      <NA>
    ## Taxa_00059                           g__Beutenbergia
    ## Taxa_00060                          g__Catenulispora
    ## Taxa_00061                           g__Cellulomonas
    ## Taxa_00062                       g__Sediminihabitans
    ## Taxa_00063                        g__Corynebacterium
    ## Taxa_00064                              g__Tomitella
    ## Taxa_00065                                      <NA>
    ## Taxa_00066                       g__Cryptosporangium
    ## Taxa_00067                             g__Fodinicola
    ## Taxa_00068                       g__Jatrophihabitans
    ## Taxa_00069                        g__Lysinimicrobium
    ## Taxa_00070                                      <NA>
    ## Taxa_00071                              g__Devriesea
    ## Taxa_00072                                      <NA>
    ## Taxa_00073                                      <NA>
    ## Taxa_00074                           g__Mobilicoccus
    ## Taxa_00075                                      <NA>
    ## Taxa_00076                           g__Blastococcus
    ## Taxa_00077                          g__Modestobacter
    ## Taxa_00078                                      <NA>
    ## Taxa_00079                         g__Arsenicicoccus
    ## Taxa_00080                            g__Phycicoccus
    ## Taxa_00081                            g__Terracoccus
    ## Taxa_00082                           g__Tetrasphaera
    ## Taxa_00083                                      <NA>
    ## Taxa_00084                          g__Angustibacter
    ## Taxa_00085                            g__Kineococcus
    ## Taxa_00086                            g__Kineosporia
    ## Taxa_00087                                      <NA>
    ## Taxa_00088                                 g__Agreia
    ## Taxa_00089                             g__Agrococcus
    ## Taxa_00090                            g__Alpinimonas
    ## Taxa_00091                          g__Amnibacterium
    ## Taxa_00092                             g__Conyzicola
    ## Taxa_00093                          g__Cryobacterium
    ## Taxa_00094                         g__Curtobacterium
    ## Taxa_00095                  g__Diaminobutyricibacter
    ## Taxa_00096                       g__Frigoribacterium
    ## Taxa_00097                         g__Frondihabitans
    ## Taxa_00098                             g__Galbitalea
    ## Taxa_00099                         g__Glaciihabitans
    ## Taxa_00100                            g__Herbiconiux
    ## Taxa_00101                       g__Homoserinibacter
    ## Taxa_00102                        g__Homoserinimonas
    ## Taxa_00103                              g__Klugiella
    ## Taxa_00104                              g__Leifsonia
    ## Taxa_00105                            g__Lysinimonas
    ## Taxa_00106                         g__Microbacterium
    ## Taxa_00107                         g__Microterricola
    ## Taxa_00108                             g__Mycetocola
    ## Taxa_00109                                 g__Naasia
    ## Taxa_00110                              g__Phycicola
    ## Taxa_00111                          g__Rathayibacter
    ## Taxa_00112                        g__Salinibacterium
    ## Taxa_00113                           g__Schumannella
    ## Taxa_00114                             g__Subtercola
    ## Taxa_00115                           g__Yonghaparkia
    ## Taxa_00116                                      <NA>
    ## Taxa_00117                             g__Acaricomes
    ## Taxa_00118                           g__Arthrobacter
    ## Taxa_00119                         g__Auritidibacter
    ## Taxa_00120                            g__Citricoccus
    ## Taxa_00121                         g__Zhihengliuella
    ## Taxa_00122                                      <NA>
    ## Taxa_00123                           g__Actinoplanes
    ## Taxa_00124                g__Allocatelliglobosispora
    ## Taxa_00125                                 g__Asanoa
    ## Taxa_00126                         g__Catellatospora
    ## Taxa_00127                    g__Catelliglobosispora
    ## Taxa_00128                          g__Couchioplanes
    ## Taxa_00129                      g__Dactylosporangium
    ## Taxa_00130                               g__Hamadaea
    ## Taxa_00131                            g__Jishengella
    ## Taxa_00132                          g__Krasilnikovia
    ## Taxa_00133                             g__Longispora
    ## Taxa_00134                          g__Luedemannella
    ## Taxa_00135                         g__Micromonospora
    ## Taxa_00136                          g__Phytohabitans
    ## Taxa_00137                              g__Pilimelia
    ## Taxa_00138                        g__Planosporangium
    ## Taxa_00139                       g__Plantactinospora
    ## Taxa_00140                        g__Polymorphospora
    ## Taxa_00141                       g__Pseudosporangium
    ## Taxa_00142                        g__Rugosimonospora
    ## Taxa_00143                         g__Spirilliplanes
    ## Taxa_00144                         g__Verrucosispora
    ## Taxa_00145                        g__Virgisporangium
    ## Taxa_00146                              g__Xiangella
    ## Taxa_00147                                      <NA>
    ## Taxa_00148                          g__Mycobacterium
    ## Taxa_00149                            g__Nakamurella
    ## Taxa_00150                                      <NA>
    ## Taxa_00151                               g__Millisia
    ## Taxa_00152                               g__Nocardia
    ## Taxa_00153                            g__Rhodococcus
    ## Taxa_00154                              g__Skermania
    ## Taxa_00155                         g__Smaragdicoccus
    ## Taxa_00156                             g__Williamsia
    ## Taxa_00157                                      <NA>
    ## Taxa_00158                          g__Aeromicrobium
    ## Taxa_00159                              g__Kribbella
    ## Taxa_00160                            g__Marmoricola
    ## Taxa_00161                                  g__Mumia
    ## Taxa_00162                           g__Nocardioides
    ## Taxa_00163                       g__Thermasporomyces
    ## Taxa_00164                                      <NA>
    ## Taxa_00165                      g__Promicromonospora
    ## Taxa_00166                        g__Xylanimicrobium
    ## Taxa_00167                                      <NA>
    ## Taxa_00168                           g__Auraticoccus
    ## Taxa_00169                         g__Friedmanniella
    ## Taxa_00170                           g__Microlunatus
    ## Taxa_00171                            g__Micropruina
    ## Taxa_00172                      g__Propionibacterium
    ## Taxa_00173                                      <NA>
    ## Taxa_00174                       g__Actinokineospora
    ## Taxa_00175                      g__Actinomycetospora
    ## Taxa_00176                        g__Actinophytocola
    ## Taxa_00177                          g__Actinosynnema
    ## Taxa_00178                      g__Alloactinosynnema
    ## Taxa_00179                          g__Amycolatopsis
    ## Taxa_00180                             g__Crossiella
    ## Taxa_00181                        g__Goodfellowiella
    ## Taxa_00182                              g__Kutzneria
    ## Taxa_00183                               g__Labedaea
    ## Taxa_00184                          g__Lechevalieria
    ## Taxa_00185                                g__Lentzea
    ## Taxa_00186                         g__Pseudonocardia
    ## Taxa_00187                          g__Saccharothrix
    ## Taxa_00188                          g__Thermobispora
    ## Taxa_00189                              g__Umezawaea
    ## Taxa_00190                            g__Yuhushiella
    ## Taxa_00191                    g__Haloactinobacterium
    ## Taxa_00192                           g__Segniliparus
    ## Taxa_00193                            g__Sporichthya
    ## Taxa_00194                                      <NA>
    ## Taxa_00195                          g__Kitasatospora
    ## Taxa_00196                      g__Streptacidiphilus
    ## Taxa_00197                           g__Streptomyces
    ## Taxa_00198                                      <NA>
    ## Taxa_00199                           g__Herbidospora
    ## Taxa_00200                             g__Nonomuraea
    ## Taxa_00201                        g__Planotetraspora
    ## Taxa_00202                      g__Sphaerisporangium
    ## Taxa_00203                      g__Streptosporangium
    ## Taxa_00204                     g__Thermocatellispora
    ## Taxa_00205                         g__Sinosporangium
    ## Taxa_00206                                      <NA>
    ## Taxa_00207                        g__Actinoallomurus
    ## Taxa_00208                         g__Actinocorallia
    ## Taxa_00209                           g__Actinomadura
    ## Taxa_00210                          g__Spirillospora
    ## Taxa_00211                        g__Thermomonospora
    ## Taxa_00212                                g__Gaiella
    ## Taxa_00213                                      <NA>
    ## Taxa_00214                           g__Conexibacter
    ## Taxa_00215                           g__Patulibacter
    ## Taxa_00216                        g__Solirubrobacter
    ## Taxa_00217                        g__Thermoleophilum
    ## Taxa_00218                                      <NA>
    ## Taxa_00219                                      <NA>
    ## Taxa_00220                    g__Thermosulfidibacter
    ## Taxa_00221                                      <NA>
    ## Taxa_00222                    g__Armatimonadetes_gp2
    ## Taxa_00223                    g__Armatimonadetes_gp4
    ## Taxa_00224                    g__Armatimonadetes_gp5
    ## Taxa_00225        g__Armatimonas/Armatimonadetes_gp1
    ## Taxa_00226       g__Chthonomonas/Armatimonadetes_gp3
    ## Taxa_00227                           g__Fimbriimonas
    ## Taxa_00228                                      <NA>
    ## Taxa_00229                                      <NA>
    ## Taxa_00230                                      <NA>
    ## Taxa_00231                                      <NA>
    ## Taxa_00232                           g__Dysgonomonas
    ## Taxa_00233                                      <NA>
    ## Taxa_00234                                      <NA>
    ## Taxa_00235                           g__Chryseolinea
    ## Taxa_00236                                      <NA>
    ## Taxa_00237                           g__Mariniradius
    ## Taxa_00238                             g__Nitritalea
    ## Taxa_00239                                      <NA>
    ## Taxa_00240                          g__Adhaeribacter
    ## Taxa_00241                              g__Cytophaga
    ## Taxa_00242                            g__Dyadobacter
    ## Taxa_00243                              g__Emticicia
    ## Taxa_00244                               g__Fibrella
    ## Taxa_00245                           g__Hymenobacter
    ## Taxa_00246                            g__Microscilla
    ## Taxa_00247                           g__Persicitalea
    ## Taxa_00248                              g__Spirosoma
    ## Taxa_00249                         g__Sporocytophaga
    ## Taxa_00250                                      <NA>
    ## Taxa_00251                            g__Aureibacter
    ## Taxa_00252                             g__Flexithrix
    ## Taxa_00253                        g__Imperialibacter
    ## Taxa_00254                             g__Limibacter
    ## Taxa_00255                          g__Marinoscillum
    ## Taxa_00256                           g__Ohtaekwangia
    ## Taxa_00257                                      <NA>
    ## Taxa_00258                                      <NA>
    ## Taxa_00259                           g__Crocinitomix
    ## Taxa_00260                             g__Fluviicola
    ## Taxa_00261                     g__Phaeocystidibacter
    ## Taxa_00262                           g__Salinirepens
    ## Taxa_00263                                      <NA>
    ## Taxa_00264                        g__Cloacibacterium
    ## Taxa_00265                         g__Flavobacterium
    ## Taxa_00266                           g__Hanstruepera
    ## Taxa_00267                             g__Imtechella
    ## Taxa_00268                            g__Namhaeicola
    ## Taxa_00269                               g__Soonwooa
    ## Taxa_00270                                      <NA>
    ## Taxa_00271                                      <NA>
    ## Taxa_00272                         g__Arachidicoccus
    ## Taxa_00273                         g__Asinibacterium
    ## Taxa_00274                           g__Chitinophaga
    ## Taxa_00275                                g__Cnuella
    ## Taxa_00276                             g__Crenotalea
    ## Taxa_00277                        g__Ferruginibacter
    ## Taxa_00278                              g__Filimonas
    ## Taxa_00279                        g__Flavihumibacter
    ## Taxa_00280                        g__Flavisolibacter
    ## Taxa_00281                             g__Flavitalea
    ## Taxa_00282                             g__Heliimonas
    ## Taxa_00283                            g__Hydrobacter
    ## Taxa_00284                             g__Lacibacter
    ## Taxa_00285                               g__Niabella
    ## Taxa_00286                              g__Niastella
    ## Taxa_00287                          g__Parafilimonas
    ## Taxa_00288                       g__Parasegetibacter
    ## Taxa_00289                      g__Sediminibacterium
    ## Taxa_00290                           g__Segetibacter
    ## Taxa_00291                             g__Taibaiella
    ## Taxa_00292                             g__Terrimonas
    ## Taxa_00293                                      <NA>
    ## Taxa_00294                      g__Haliscomenobacter
    ## Taxa_00295                     g__Phaeodactylibacter
    ## Taxa_00296                            g__Portibacter
    ## Taxa_00297                                      <NA>
    ## Taxa_00298                          g__Arcticibacter
    ## Taxa_00299                       g__Mucilaginibacter
    ## Taxa_00300                               g__Nubsella
    ## Taxa_00301                             g__Pedobacter
    ## Taxa_00302                 g__Pseudosphingobacterium
    ## Taxa_00303                              g__Solitalea
    ## Taxa_00304             g__BRC1_genera_incertae_sedis
    ## Taxa_00305            g__WPS-1_genera_incertae_sedis
    ## Taxa_00306            g__WPS-2_genera_incertae_sedis
    ## Taxa_00307              g__ZB3_genera_incertae_sedis
    ## Taxa_00308 g__Saccharibacteria_genera_incertae_sedis
    ## Taxa_00309                                      <NA>
    ## Taxa_00310                                      <NA>
    ## Taxa_00311                           g__Neochlamydia
    ## Taxa_00312                          g__Parachlamydia
    ## Taxa_00313                               g__Simkania
    ## Taxa_00314                                      <NA>
    ## Taxa_00315                                      <NA>
    ## Taxa_00316                             g__Bellilinea
    ## Taxa_00317                             g__Leptolinea
    ## Taxa_00318                              g__Levilinea
    ## Taxa_00319                             g__Longilinea
    ## Taxa_00320                              g__Pelolinea
    ## Taxa_00321                          g__Ardenticatena
    ## Taxa_00322                                      <NA>
    ## Taxa_00323                             g__Caldilinea
    ## Taxa_00324                            g__Litorilinea
    ## Taxa_00325                                      <NA>
    ## Taxa_00326                                      <NA>
    ## Taxa_00327                                      <NA>
    ## Taxa_00328                             g__Heliothrix
    ## Taxa_00329                            g__Roseiflexus
    ## Taxa_00330                             g__Kallotenue
    ## Taxa_00331                        g__Dehalococcoides
    ## Taxa_00332                                      <NA>
    ## Taxa_00333                                      <NA>
    ## Taxa_00334                          g__Ktedonobacter
    ## Taxa_00335                       g__Thermosporothrix
    ## Taxa_00336                           g__Thermoflexus
    ## Taxa_00337                                      <NA>
    ## Taxa_00338                                      <NA>
    ## Taxa_00339                          g__Sphaerobacter
    ## Taxa_00340                            g__Thermorudis
    ## Taxa_00341                                      <NA>
    ## Taxa_00342                                      <NA>
    ## Taxa_00343                        g__Bacillariophyta
    ## Taxa_00344                          g__Bangiophyceae
    ## Taxa_00345                   g__Chlorarachniophyceae
    ## Taxa_00346                            g__Chlorophyta
    ## Taxa_00347                           g__Streptophyta
    ## Taxa_00348                                    g__GpI
    ## Taxa_00349                                      <NA>
    ## Taxa_00350                                      <NA>
    ## Taxa_00351                         g__Elusimicrobium
    ## Taxa_00352               g__Candidatus Endomicrobium
    ## Taxa_00353                                      <NA>
    ## Taxa_00354                                      <NA>
    ## Taxa_00355                                      <NA>
    ## Taxa_00356                                      <NA>
    ## Taxa_00357                       g__Alicyclobacillus
    ## Taxa_00358                         g__Effusibacillus
    ## Taxa_00359                               g__Kyrpidia
    ## Taxa_00360                           g__Tumebacillus
    ## Taxa_00361                                      <NA>
    ## Taxa_00362                           g__Aeribacillus
    ## Taxa_00363                               g__Bacillus
    ## Taxa_00364                           g__Domibacillus
    ## Taxa_00365                          g__Falsibacillus
    ## Taxa_00366                            g__Geobacillus
    ## Taxa_00367                         g__Saccharococcus
    ## Taxa_00368                                      <NA>
    ## Taxa_00369                           g__Allobacillus
    ## Taxa_00370                         g__Cerasibacillus
    ## Taxa_00371                       g__Compostibacillus
    ## Taxa_00372                      g__Halalkalibacillus
    ## Taxa_00373                       g__Melghiribacillus
    ## Taxa_00374                        g__Natronobacillus
    ## Taxa_00375                         g__Oceanobacillus
    ## Taxa_00376                      g__Paucisalibacillus
    ## Taxa_00377                      g__Pullulanibacillus
    ## Taxa_00378                            g__Salirhabdus
    ## Taxa_00379                      g__Saliterribacillus
    ## Taxa_00380                         g__Tuberibacillus
    ## Taxa_00381                            g__Thermicanus
    ## Taxa_00382                                      <NA>
    ## Taxa_00383                      g__Desulfuribacillus
    ## Taxa_00384                                      <NA>
    ## Taxa_00385                          g__Brevibacillus
    ## Taxa_00386                               g__Cohnella
    ## Taxa_00387                          g__Fontibacillus
    ## Taxa_00388                          g__Paenibacillus
    ## Taxa_00389                       g__Saccharibacillus
    ## Taxa_00390                         g__Thermobacillus
    ## Taxa_00391                                      <NA>
    ## Taxa_00392                           g__Ammoniphilus
    ## Taxa_00393                       g__Aneurinibacillus
    ## Taxa_00394                            g__Oxalophagus
    ## Taxa_00395                                      <NA>
    ## Taxa_00396                            g__Caryophanon
    ## Taxa_00397                             g__Chungangia
    ## Taxa_00398                             g__Filibacter
    ## Taxa_00399                       g__Jeotgalibacillus
    ## Taxa_00400                         g__Lysinibacillus
    ## Taxa_00401                      g__Paenisporosarcina
    ## Taxa_00402          g__Planococcaceae_incertae_sedis
    ## Taxa_00403                        g__Psychrobacillus
    ## Taxa_00404                       g__Rummeliibacillus
    ## Taxa_00405                           g__Solibacillus
    ## Taxa_00406                           g__Sporosarcina
    ## Taxa_00407                         g__Viridibacillus
    ## Taxa_00408                         g__Staphylococcus
    ## Taxa_00409                                      <NA>
    ## Taxa_00410                             g__Desmospora
    ## Taxa_00411                            g__Lihuaxuella
    ## Taxa_00412                            g__Shimazuella
    ## Taxa_00413                      g__Thermoactinomyces
    ## Taxa_00414                   g__Thermoflavimicrobium
    ## Taxa_00415                                      <NA>
    ## Taxa_00416                             g__Planifilum
    ## Taxa_00417                         g__Polycladomyces
    ## Taxa_00418                                      <NA>
    ## Taxa_00419                               g__Desemzia
    ## Taxa_00420                             g__Pilibacter
    ## Taxa_00421                          g__Lactobacillus
    ## Taxa_00422                            g__Lactococcus
    ## Taxa_00423                          g__Streptococcus
    ## Taxa_00424                                      <NA>
    ## Taxa_00425                                      <NA>
    ## Taxa_00426                             g__Catabacter
    ## Taxa_00427                                      <NA>
    ## Taxa_00428                           g__Anaerobacter
    ## Taxa_00429              g__Clostridium sensu stricto
    ## Taxa_00430                           g__Fervidicella
    ## Taxa_00431                          g__Oceanirhabdus
    ## Taxa_00432                              g__Oxobacter
    ## Taxa_00433                      g__Proteiniclasticum
    ## Taxa_00434                                      <NA>
    ## Taxa_00435                      g__Tepidanaerobacter
    ## Taxa_00436                  g__Thermoanaerobacterium
    ## Taxa_00437                           g__Anaerococcus
    ## Taxa_00438                            g__Anaerovorax
    ## Taxa_00439                        g__Symbiobacterium
    ## Taxa_00440                                      <NA>
    ## Taxa_00441                          g__Gracilibacter
    ## Taxa_00442                         g__Hydrogenispora
    ## Taxa_00443                                      <NA>
    ## Taxa_00444                       g__Clostridium XlVa
    ## Taxa_00445                         g__Eisenbergiella
    ## Taxa_00446                            g__Mobilitalea
    ## Taxa_00447                                      <NA>
    ## Taxa_00448                      g__Desulfosporosinus
    ## Taxa_00449                                      <NA>
    ## Taxa_00450                    g__Acetanaerobacterium
    ## Taxa_00451                            g__Acetivibrio
    ## Taxa_00452                        g__Anaerobacterium
    ## Taxa_00453                        g__Cellulosibacter
    ## Taxa_00454                        g__Clostridium III
    ## Taxa_00455                         g__Clostridium IV
    ## Taxa_00456                         g__Ethanoligenens
    ## Taxa_00457                      g__Pseudobacteroides
    ## Taxa_00458                            g__Sporobacter
    ## Taxa_00459                      g__Thermohydrogenium
    ## Taxa_00460                      g__Thermodesulfobium
    ## Taxa_00461                                      <NA>
    ## Taxa_00462                                      <NA>
    ## Taxa_00463                           g__Succinispira
    ## Taxa_00464                                      <NA>
    ## Taxa_00465                              g__Pelosinus
    ## Taxa_00466                           g__Psychrosinus
    ## Taxa_00467                           g__Gemmatimonas
    ## Taxa_00468               g__Candidatus Hydrogenedens
    ## Taxa_00469  g__Latescibacteria_genera_incertae_sedis
    ## Taxa_00470                                      <NA>
    ## Taxa_00471                            g__Victivallis
    ## Taxa_00472   g__Microgenomates_genera_incertae_sedis
    ## Taxa_00473                             g__Nitrospira
    ## Taxa_00474     g__Omnitrophica_genera_incertae_sedis
    ## Taxa_00475    g__Parcubacteria_genera_incertae_sedis
    ## Taxa_00476                                      <NA>
    ## Taxa_00477                                      <NA>
    ## Taxa_00478                           g__Phycisphaera
    ## Taxa_00479                          g__Tepidisphaera
    ## Taxa_00480                                      <NA>
    ## Taxa_00481                                      <NA>
    ## Taxa_00482              g__Candidatus Anammoxoglobus
    ## Taxa_00483                                      <NA>
    ## Taxa_00484                            g__Aquisphaera
    ## Taxa_00485                        g__Blastopirellula
    ## Taxa_00486                                g__Gemmata
    ## Taxa_00487                                g__Gimesia
    ## Taxa_00488                             g__Isosphaera
    ## Taxa_00489                              g__Pirellula
    ## Taxa_00490                       g__Planctomicrobium
    ## Taxa_00491                           g__Planctopirus
    ## Taxa_00492                          g__Rubinisphaera
    ## Taxa_00493                            g__Schlesneria
    ## Taxa_00494                         g__Singulisphaera
    ## Taxa_00495                            g__Telmatocola
    ## Taxa_00496                            g__Thermogutta
    ## Taxa_00497                           g__Zavarzinella
    ## Taxa_00498                                      <NA>
    ## Taxa_00499                                      <NA>
    ## Taxa_00500                                      <NA>
    ## Taxa_00501                             g__Breoghania
    ## Taxa_00502                           g__Geminicoccus
    ## Taxa_00503                         g__Rhizomicrobium
    ## Taxa_00504                                      <NA>
    ## Taxa_00505                                      <NA>
    ## Taxa_00506                          g__Asticcacaulis
    ## Taxa_00507                          g__Brevundimonas
    ## Taxa_00508                            g__Caulobacter
    ## Taxa_00509                       g__Phenylobacterium
    ## Taxa_00510                                      <NA>
    ## Taxa_00511                             g__Maricaulis
    ## Taxa_00512                            g__Marinicauda
    ## Taxa_00513                            g__Eilatimonas
    ## Taxa_00514                          g__Amphiplicatus
    ## Taxa_00515                                      <NA>
    ## Taxa_00516                             g__Aureimonas
    ## Taxa_00517                                      <NA>
    ## Taxa_00518                           g__Beijerinckia
    ## Taxa_00519                            g__Camelimonas
    ## Taxa_00520                          g__Chelatococcus
    ## Taxa_00521                           g__Methylocapsa
    ## Taxa_00522                           g__Methylocella
    ## Taxa_00523                          g__Methyloferula
    ## Taxa_00524                          g__Methylorosula
    ## Taxa_00525                         g__Methylovirgula
    ## Taxa_00526                    g__Pseudochelatococcus
    ## Taxa_00527                                      <NA>
    ## Taxa_00528                                 g__Afipia
    ## Taxa_00529                                  g__Bosea
    ## Taxa_00530                         g__Bradyrhizobium
    ## Taxa_00531                            g__Nitrobacter
    ## Taxa_00532                            g__Oligotropha
    ## Taxa_00533                           g__Rhodoblastus
    ## Taxa_00534                       g__Rhodopseudomonas
    ## Taxa_00535                          g__Salinarimonas
    ## Taxa_00536                             g__Tardiphaga
    ## Taxa_00537                                      <NA>
    ## Taxa_00538                                g__Daeguia
    ## Taxa_00539                       g__Falsochrobactrum
    ## Taxa_00540                       g__Paenochrobactrum
    ## Taxa_00541                                      <NA>
    ## Taxa_00542                        g__Ancalomicrobium
    ## Taxa_00543                             g__Aquabacter
    ## Taxa_00544                          g__Blastochloris
    ## Taxa_00545                           g__Cucumibacter
    ## Taxa_00546                                g__Devosia
    ## Taxa_00547                          g__Filomicrobium
    ## Taxa_00548                         g__Hyphomicrobium
    ## Taxa_00549                          g__Pedomicrobium
    ## Taxa_00550                     g__Prosthecomicrobium
    ## Taxa_00551                         g__Rhodomicrobium
    ## Taxa_00552                            g__Rhodoplanes
    ## Taxa_00553                                      <NA>
    ## Taxa_00554                               g__Meganema
    ## Taxa_00555                       g__Methylobacterium
    ## Taxa_00556                             g__Microvirga
    ## Taxa_00557                      g__Psychroglaciecola
    ## Taxa_00558                                      <NA>
    ## Taxa_00559                             g__Albibacter
    ## Taxa_00560                          g__Methylocystis
    ## Taxa_00561                            g__Methylopila
    ## Taxa_00562                           g__Methylosinus
    ## Taxa_00563                        g__Pleomorphomonas
    ## Taxa_00564                                      <NA>
    ## Taxa_00565                            g__Aminobacter
    ## Taxa_00566                          g__Aquamicrobium
    ## Taxa_00567                          g__Chelativorans
    ## Taxa_00568                                g__Hoeflea
    ## Taxa_00569                      g__Lentilitoribacter
    ## Taxa_00570                          g__Mesorhizobium
    ## Taxa_00571                        g__Phyllobacterium
    ## Taxa_00572                       g__Pseudaminobacter
    ## Taxa_00573                                      <NA>
    ## Taxa_00574                           g__Ciceribacter
    ## Taxa_00575                                g__Ensifer
    ## Taxa_00576                                g__Kaistia
    ## Taxa_00577                           g__Neorhizobium
    ## Taxa_00578                              g__Rhizobium
    ## Taxa_00579                               g__Shinella
    ## Taxa_00580                                      <NA>
    ## Taxa_00581                             g__Alsobacter
    ## Taxa_00582                                g__Bauldia
    ## Taxa_00583                          g__Phreatobacter
    ## Taxa_00584                            g__Variibacter
    ## Taxa_00585                            g__Vasilyevaea
    ## Taxa_00586                                      <NA>
    ## Taxa_00587                               g__Afifella
    ## Taxa_00588                       g__Dichotomicrobium
    ## Taxa_00589                            g__Lutibaculum
    ## Taxa_00590                     g__Methyloceanibacter
    ## Taxa_00591                         g__Methyloligella
    ## Taxa_00592                           g__Parvibaculum
    ## Taxa_00593                       g__Rhodoligotrophos
    ## Taxa_00594                          g__Tepidamorphus
    ## Taxa_00595                             g__Roseiarcus
    ## Taxa_00596                                      <NA>
    ## Taxa_00597                                 g__Labrys
    ## Taxa_00598                           g__Pseudolabrys
    ## Taxa_00599                                      <NA>
    ## Taxa_00600                            g__Agaricicola
    ## Taxa_00601                               g__Ahrensia
    ## Taxa_00602                            g__Cereibacter
    ## Taxa_00603                          g__Frigidibacter
    ## Taxa_00604                            g__Gemmobacter
    ## Taxa_00605                           g__Hasllibacter
    ## Taxa_00606                      g__Pseudorhodobacter
    ## Taxa_00607                                      <NA>
    ## Taxa_00608                                      <NA>
    ## Taxa_00609                            g__Acetobacter
    ## Taxa_00610                              g__Acidisoma
    ## Taxa_00611                           g__Acidisphaera
    ## Taxa_00612                             g__Acidocella
    ## Taxa_00613                             g__Acidomonas
    ## Taxa_00614                              g__Ameyamaea
    ## Taxa_00615                                  g__Asaia
    ## Taxa_00616                           g__Craurococcus
    ## Taxa_00617                             g__Endobacter
    ## Taxa_00618                      g__Gluconacetobacter
    ## Taxa_00619                          g__Granulibacter
    ## Taxa_00620                              g__Humitalea
    ## Taxa_00621                                g__Kozakia
    ## Taxa_00622                               g__Neoasaia
    ## Taxa_00623                          g__Neokomagataea
    ## Taxa_00624                       g__Paracraurococcus
    ## Taxa_00625                              g__Rhodopila
    ## Taxa_00626                            g__Rhodovarius
    ## Taxa_00627                            g__Roseococcus
    ## Taxa_00628                             g__Roseomonas
    ## Taxa_00629                                 g__Stella
    ## Taxa_00630                          g__Swaminathania
    ## Taxa_00631                         g__Tanticharoenia
    ## Taxa_00632                               g__Elioraea
    ## Taxa_00633                             g__Reyranella
    ## Taxa_00634                                      <NA>
    ## Taxa_00635                           g__Azospirillum
    ## Taxa_00636                       g__Constrictibacter
    ## Taxa_00637                         g__Defluviicoccus
    ## Taxa_00638                                 g__Dongia
    ## Taxa_00639                             g__Inquilinus
    ## Taxa_00640                          g__Lacibacterium
    ## Taxa_00641                              g__Limimonas
    ## Taxa_00642                       g__Magnetospirillum
    ## Taxa_00643                          g__Marispirillum
    ## Taxa_00644                                 g__Nisaea
    ## Taxa_00645                         g__Nitrospirillum
    ## Taxa_00646                         g__Niveispirillum
    ## Taxa_00647                          g__Novispirillum
    ## Taxa_00648                          g__Oceanibaculum
    ## Taxa_00649                         g__Phaeospirillum
    ## Taxa_00650                            g__Skermanella
    ## Taxa_00651                       g__Telmatospirillum
    ## Taxa_00652                                g__Tistlia
    ## Taxa_00653                                      <NA>
    ## Taxa_00654                              g__Anaplasma
    ## Taxa_00655                                      <NA>
    ## Taxa_00656                               g__Orientia
    ## Taxa_00657                             g__Rickettsia
    ## Taxa_00658                                      <NA>
    ## Taxa_00659                               g__Taonella
    ## Taxa_00660                                      <NA>
    ## Taxa_00661                                      <NA>
    ## Taxa_00662                     g__Altererythrobacter
    ## Taxa_00663                       g__Erythromicrobium
    ## Taxa_00664                         g__Porphyrobacter
    ## Taxa_00665                                      <NA>
    ## Taxa_00666                            g__Blastomonas
    ## Taxa_00667                             g__Hephaestia
    ## Taxa_00668                        g__Novosphingobium
    ## Taxa_00669                        g__Parablastomonas
    ## Taxa_00670                       g__Parasphingopyxis
    ## Taxa_00671                       g__Polymorphobacter
    ## Taxa_00672                           g__Rhizorhabdus
    ## Taxa_00673                            g__Rhizorhapis
    ## Taxa_00674                      g__Sandaracinobacter
    ## Taxa_00675                     g__Sandarakinorhabdus
    ## Taxa_00676                            g__Sphingobium
    ## Taxa_00677                       g__Sphingomicrobium
    ## Taxa_00678                           g__Sphingomonas
    ## Taxa_00679                           g__Sphingopyxis
    ## Taxa_00680                         g__Sphingorhabdus
    ## Taxa_00681                       g__Sphingosinicella
    ## Taxa_00682                                      <NA>
    ## Taxa_00683                                      <NA>
    ## Taxa_00684                                      <NA>
    ## Taxa_00685                          g__Achromobacter
    ## Taxa_00686                            g__Alcaligenes
    ## Taxa_00687                           g__Candidimonas
    ## Taxa_00688                                 g__Derxia
    ## Taxa_00689                             g__Eoetvoesia
    ## Taxa_00690                        g__Paenalcaligenes
    ## Taxa_00691                         g__Paralcaligenes
    ## Taxa_00692                          g__Pigmentiphaga
    ## Taxa_00693                           g__Pusillimonas
    ## Taxa_00694                                      <NA>
    ## Taxa_00695                           g__Burkholderia
    ## Taxa_00696                           g__Chitinimonas
    ## Taxa_00697                              g__Lautropia
    ## Taxa_00698                            g__Limnobacter
    ## Taxa_00699                             g__Paucimonas
    ## Taxa_00700                                      <NA>
    ## Taxa_00701                          g__Aquabacterium
    ## Taxa_00702                              g__Aquincola
    ## Taxa_00703                              g__Ideonella
    ## Taxa_00704                             g__Leptothrix
    ## Taxa_00705                            g__Methylibium
    ## Taxa_00706                            g__Paucibacter
    ## Taxa_00707                          g__Piscinibacter
    ## Taxa_00708                             g__Rivibacter
    ## Taxa_00709                             g__Rubrivivax
    ## Taxa_00710                             g__Thiobacter
    ## Taxa_00711                             g__Xylophilus
    ## Taxa_00712                                      <NA>
    ## Taxa_00713                             g__Acidovorax
    ## Taxa_00714                         g__Alicycliphilus
    ## Taxa_00715                             g__Caenimonas
    ## Taxa_00716                              g__Comamonas
    ## Taxa_00717                            g__Curvibacter
    ## Taxa_00718                         g__Diaphorobacter
    ## Taxa_00719                         g__Hydrogenophaga
    ## Taxa_00720                            g__Hylemonella
    ## Taxa_00721                             g__Kinneretia
    ## Taxa_00722                                g__Malikia
    ## Taxa_00723                                g__Ottowia
    ## Taxa_00724                              g__Pelomonas
    ## Taxa_00725                            g__Polaromonas
    ## Taxa_00726                        g__Pseudacidovorax
    ## Taxa_00727                       g__Pseudorhodoferax
    ## Taxa_00728                            g__Ramlibacter
    ## Taxa_00729                             g__Rhodoferax
    ## Taxa_00730                           g__Schlegelella
    ## Taxa_00731                          g__Simplicispira
    ## Taxa_00732                             g__Variovorax
    ## Taxa_00733                             g__Xenophilus
    ## Taxa_00734                         g__Zhizhongheella
    ## Taxa_00735                                      <NA>
    ## Taxa_00736                             g__Collimonas
    ## Taxa_00737                              g__Duganella
    ## Taxa_00738                            g__Glaciimonas
    ## Taxa_00739                         g__Herbaspirillum
    ## Taxa_00740                          g__Herminiimonas
    ## Taxa_00741                      g__Janthinobacterium
    ## Taxa_00742                               g__Massilia
    ## Taxa_00743                     g__Noviherbaspirillum
    ## Taxa_00744                       g__Oxalicibacterium
    ## Taxa_00745                     g__Paraherbaspirillum
    ## Taxa_00746                        g__Pseudoduganella
    ## Taxa_00747                              g__Rugamonas
    ## Taxa_00748                          g__Undibacterium
    ## Taxa_00749                             g__Sutterella
    ## Taxa_00750                         g__Ferritrophicum
    ## Taxa_00751                               g__Ferrovum
    ## Taxa_00752                                      <NA>
    ## Taxa_00753                            g__Gallionella
    ## Taxa_00754                           g__Sideroxydans
    ## Taxa_00755                                      <NA>
    ## Taxa_00756                          g__Methylotenera
    ## Taxa_00757                                      <NA>
    ## Taxa_00758                         g__Amantichitinum
    ## Taxa_00759                          g__Aquaspirillum
    ## Taxa_00760                          g__Chitiniphilus
    ## Taxa_00761                                  g__Leeia
    ## Taxa_00762                            g__Simonsiella
    ## Taxa_00763                          g__Snodgrassella
    ## Taxa_00764                                      <NA>
    ## Taxa_00765                           g__Nitrosospira
    ## Taxa_00766                 g__Candidatus Procabacter
    ## Taxa_00767                                      <NA>
    ## Taxa_00768                              g__Azovibrio
    ## Taxa_00769                          g__Denitratisoma
    ## Taxa_00770                           g__Georgfuchsia
    ## Taxa_00771                      g__Methyloversatilis
    ## Taxa_00772                            g__Rhodocyclus
    ## Taxa_00773                       g__Sterolibacterium
    ## Taxa_00774                            g__Sulfurisoma
    ## Taxa_00775                           g__Sulfuritalea
    ## Taxa_00776                     g__Uliginosibacterium
    ## Taxa_00777                                      <NA>
    ## Taxa_00778                           g__Sulfuricella
    ## Taxa_00779                                      <NA>
    ## Taxa_00780                                      <NA>
    ## Taxa_00781                                      <NA>
    ## Taxa_00782                          g__Bacteriovorax
    ## Taxa_00783                           g__Peredibacter
    ## Taxa_00784                                      <NA>
    ## Taxa_00785                           g__Bdellovibrio
    ## Taxa_00786                          g__Vampirovibrio
    ## Taxa_00787                    g__Pseudobacteriovorax
    ## Taxa_00788                                      <NA>
    ## Taxa_00789                            g__Deferrisoma
    ## Taxa_00790                       g__Dissulfuribacter
    ## Taxa_00791                                      <NA>
    ## Taxa_00792                                      <NA>
    ## Taxa_00793                         g__Desulfatitalea
    ## Taxa_00794                                      <NA>
    ## Taxa_00795                         g__Desulfobaculum
    ## Taxa_00796                                      <NA>
    ## Taxa_00797                                      <NA>
    ## Taxa_00798                              g__Geobacter
    ## Taxa_00799                       g__Geopsychrobacter
    ## Taxa_00800                                      <NA>
    ## Taxa_00801                                      <NA>
    ## Taxa_00802                       g__Anaeromyxobacter
    ## Taxa_00803                             g__Archangium
    ## Taxa_00804                            g__Cystobacter
    ## Taxa_00805                             g__Hyalangium
    ## Taxa_00806                           g__Melittangium
    ## Taxa_00807                            g__Stigmatella
    ## Taxa_00808                             g__Haliangium
    ## Taxa_00809                               g__Kofleria
    ## Taxa_00810                            g__Labilithrix
    ## Taxa_00811                                      <NA>
    ## Taxa_00812                          g__Aggregicoccus
    ## Taxa_00813                             g__Myxococcus
    ## Taxa_00814                           g__Pyxidicoccus
    ## Taxa_00815                                      <NA>
    ## Taxa_00816                            g__Nannocystis
    ## Taxa_00817                           g__Plesiocystis
    ## Taxa_00818                          g__Phaselicystis
    ## Taxa_00819                                      <NA>
    ## Taxa_00820                             g__Byssovorax
    ## Taxa_00821                               g__Jahnella
    ## Taxa_00822                             g__Minicystis
    ## Taxa_00823                             g__Polyangium
    ## Taxa_00824                              g__Sorangium
    ## Taxa_00825                           g__Sandaracinus
    ## Taxa_00826                          g__Vulgatibacter
    ## Taxa_00827                                      <NA>
    ## Taxa_00828                              g__Smithella
    ## Taxa_00829                                      <NA>
    ## Taxa_00830                                      <NA>
    ## Taxa_00831                                  g__Cetia
    ## Taxa_00832                                      <NA>
    ## Taxa_00833                                      <NA>
    ## Taxa_00834                                      <NA>
    ## Taxa_00835                             g__Chromatium
    ## Taxa_00836                          g__Nitrosococcus
    ## Taxa_00837                               g__Thiobaca
    ## Taxa_00838                        g__Thioflavicoccus
    ## Taxa_00839                         g__Thiohalobacter
    ## Taxa_00840                                      <NA>
    ## Taxa_00841                       g__Acidiferrobacter
    ## Taxa_00842                     g__Ectothiorhodosinus
    ## Taxa_00843                           g__Natronocella
    ## Taxa_00844                            g__Spiribacter
    ## Taxa_00845                              g__Thioalbus
    ## Taxa_00846                        g__Thioalkalispira
    ## Taxa_00847                         g__Thiorhodospira
    ## Taxa_00848                                      <NA>
    ## Taxa_00849                              g__Ewingella
    ## Taxa_00850                        g__Obesumbacterium
    ## Taxa_00851                                g__Pantoea
    ## Taxa_00852                               g__Yersinia
    ## Taxa_00853                                      <NA>
    ## Taxa_00854                  g__Candidatus Carsonella
    ## Taxa_00855                       g__Methylohalomonas
    ## Taxa_00856                          g__Sedimenticola
    ## Taxa_00857                                      <NA>
    ## Taxa_00858                                      <NA>
    ## Taxa_00859                              g__Aquicella
    ## Taxa_00860                               g__Coxiella
    ## Taxa_00861                        g__Diplorickettsia
    ## Taxa_00862                             g__Legionella
    ## Taxa_00863                                      <NA>
    ## Taxa_00864                          g__Methylococcus
    ## Taxa_00865                            g__Methylogaea
    ## Taxa_00866                         g__Methylomarinum
    ## Taxa_00867                      g__Methyloparacoccus
    ## Taxa_00868                                      <NA>
    ## Taxa_00869                            g__Zooshikella
    ## Taxa_00870                          g__Motiliproteus
    ## Taxa_00871                        g__Pseudohongiella
    ## Taxa_00872                            g__Haemophilus
    ## Taxa_00873                                      <NA>
    ## Taxa_00874                           g__Alkanindiges
    ## Taxa_00875                                      <NA>
    ## Taxa_00876                               g__Azomonas
    ## Taxa_00877                             g__Cellvibrio
    ## Taxa_00878                            g__Pseudomonas
    ## Taxa_00879                            g__Rhizobacter
    ## Taxa_00880                                g__Serpens
    ## Taxa_00881                                      <NA>
    ## Taxa_00882                              g__Beggiatoa
    ## Taxa_00883                                      <NA>
    ## Taxa_00884                             g__Algiphilus
    ## Taxa_00885                                      <NA>
    ## Taxa_00886                             g__Fontimonas
    ## Taxa_00887                      g__Hydrocarboniphaga
    ## Taxa_00888                                g__Nevskia
    ## Taxa_00889                         g__Panacagrimonas
    ## Taxa_00890                           g__Povalibacter
    ## Taxa_00891                              g__Solimonas
    ## Taxa_00892                         g__Steroidobacter
    ## Taxa_00893                                      <NA>
    ## Taxa_00894                              g__Aquimonas
    ## Taxa_00895                             g__Arenimonas
    ## Taxa_00896                             g__Aspromonas
    ## Taxa_00897                           g__Chiayiivirga
    ## Taxa_00898                             g__Dokdonella
    ## Taxa_00899                                 g__Dyella
    ## Taxa_00900                              g__Frateuria
    ## Taxa_00901                            g__Luteibacter
    ## Taxa_00902                             g__Lysobacter
    ## Taxa_00903                       g__Metallibacterium
    ## Taxa_00904                        g__Mizugakiibacter
    ## Taxa_00905                      g__Pseudoxanthomonas
    ## Taxa_00906                         g__Rehaibacterium
    ## Taxa_00907                          g__Rhodanobacter
    ## Taxa_00908                                 g__Rudaea
    ## Taxa_00909                       g__Stenotrophomonas
    ## Taxa_00910                             g__Tahibacter
    ## Taxa_00911                      g__Vulcaniibacterium
    ## Taxa_00912                            g__Xanthomonas
    ## Taxa_00913                                g__Xylella
    ## Taxa_00914                            g__Oligoflexus
    ## Taxa_00915                                      <NA>
    ## Taxa_00916                              g__Brevinema
    ## Taxa_00917                                      <NA>
    ## Taxa_00918                              g__Leptonema
    ## Taxa_00919                            g__Turneriella
    ## Taxa_00920                                      <NA>
    ## Taxa_00921                            g__Salinispira
    ## Taxa_00922                            g__Spirochaeta
    ## Taxa_00923              g__SR1_genera_incertae_sedis
    ## Taxa_00924                           g__Entomoplasma
    ## Taxa_00925                             g__Haloplasma
    ## Taxa_00926                                      <NA>
    ## Taxa_00927                                      <NA>
    ## Taxa_00928                                      <NA>
    ## Taxa_00929                                      <NA>
    ## Taxa_00930                           g__Alterococcus
    ## Taxa_00931                               g__Opitutus
    ## Taxa_00932                                      <NA>
    ## Taxa_00933                           g__Cerasicoccus
    ## Taxa_00934                          g__Puniceicoccus
    ## Taxa_00935                                      <NA>
    ## Taxa_00936   g__Spartobacteria_genera_incertae_sedis
    ## Taxa_00937                         g__Terrimicrobium
    ## Taxa_00938                                      <NA>
    ## Taxa_00939                            g__Limisphaera
    ## Taxa_00940     g__Subdivision3_genera_incertae_sedis
    ## Taxa_00941                                      <NA>
    ## Taxa_00942                            g__Brevifollis
    ## Taxa_00943                             g__Haloferula
    ## Taxa_00944                          g__Luteolibacter
    ## Taxa_00945                         g__Persicirhabdus
    ## Taxa_00946                        g__Prosthecobacter
    ## Taxa_00947                         g__Roseimicrobium
    ## Taxa_00948                       g__Verrucomicrobium
    ## Taxa_00949                                      <NA>

The following ensures that features with ambiguous phylum annotation are also removed.

``` r
s16sV1V3.1 <- subset_taxa(s16sV1V3, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
s16sV1V3.1
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 947 taxa and 28 samples ]
    ## sample_data() Sample Data:       [ 28 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 947 taxa by 6 taxonomic ranks ]

### Now subset low occurance taxa

``` r
plyr::ddply(prevelancedf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),mean_abundance=mean(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
  })
```

    ##                            Phylum mean_prevalence mean_abundance
    ## 1                p__Acidobacteria       20.463415    3494.317073
    ## 2               p__Actinobacteria       10.103448     298.804598
    ## 3                    p__Aquificae        1.333333       1.333333
    ## 4              p__Armatimonadetes       23.857143     623.285714
    ## 5                p__Bacteroidetes       12.171053     471.144737
    ## 6                         p__BRC1       26.000000      95.000000
    ## 7     p__candidate division WPS-1       28.000000    6988.000000
    ## 8     p__candidate division WPS-2       27.000000    1083.000000
    ## 9       p__candidate division ZB3        1.000000       1.000000
    ## 10 p__Candidatus Saccharibacteria       28.000000    3335.000000
    ## 11                  p__Chlamydiae        8.800000      16.200000
    ## 12                 p__Chloroflexi       13.444444     254.851852
    ## 13               p__Crenarchaeota        1.000000       1.000000
    ## 14   p__Cyanobacteria/Chloroplast        8.625000      76.625000
    ## 15         p__Deinococcus-Thermus        3.000000       4.000000
    ## 16               p__Elusimicrobia        9.333333      16.666667
    ## 17                  p__Firmicutes        8.798246     519.017544
    ## 18            p__Gemmatimonadetes       28.000000   22581.000000
    ## 19             p__Hydrogenedentes        9.000000      16.000000
    ## 20             p__Latescibacteria       25.000000     398.000000
    ## 21               p__Lentisphaerae        5.000000      13.500000
    ## 22              p__Microgenomates       21.000000     192.000000
    ## 23                 p__Nitrospirae       22.000000    1387.000000
    ## 24                p__Omnitrophica        3.000000       3.000000
    ## 25               p__Parcubacteria       28.000000    5848.000000
    ## 26              p__Planctomycetes       21.363636    1149.727273
    ## 27              p__Proteobacteria       12.359712     728.748201
    ## 28                p__Spirochaetes        5.125000      10.625000
    ## 29                         p__SR1        7.000000      15.000000
    ## 30                 p__Tenericutes        1.000000       1.000000
    ## 31       p__Thermodesulfobacteria        1.000000       1.000000
    ## 32             p__Verrucomicrobia       19.272727    2488.545455
    ## 33                           <NA>       21.333333   21367.333333

### Define phyla to filter

``` r
phyla2Filter = c("p__Aquificae", "p__candidate division ZB3",
  "p__Crenarchaeota","p__Deinococcus-Thermus","p__Omnitrophica","p__Tenericutes","p__Thermodesulfobacteria")
# Filter entries with unidentified Phylum.
s16sV1V3.2 = subset_taxa(s16sV1V3.1, !Phylum %in% phyla2Filter)
s16sV1V3.2
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 937 taxa and 28 samples ]
    ## sample_data() Sample Data:       [ 28 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 937 taxa by 6 taxonomic ranks ]

### Subset to the remaining phyla

``` r
prevelancedf1 = subset(prevelancedf, Phylum %in% get_taxa_unique(s16sV1V3.2, taxonomic.rank = "Phylum"))
ggplot(prevelancedf1, aes(TotalAbundance, Prevalence / nsamples(s16sV1V3.2),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
```

![](phyloseq_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
#  Define prevalence threshold as 10% of total samples
prevalenceThreshold = 0.10 * nsamples(s16sV1V3.2)
prevalenceThreshold
```

    ## [1] 2.8

``` r
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevelancedf1)[(prevelancedf1$Prevalence >= prevalenceThreshold)]
s16sV1V3.3 = prune_taxa(keepTaxa, s16sV1V3.2)
s16sV1V3.3
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 710 taxa and 28 samples ]
    ## sample_data() Sample Data:       [ 28 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 710 taxa by 6 taxonomic ranks ]

### Agglomerate taxa and remove all taxa with not at the genus level

``` r
length(get_taxa_unique(s16sV1V3.3, taxonomic.rank = "Genus"))
```

    ## [1] 552

``` r
s16sV1V3.4 = tax_glom(s16sV1V3.3, "Genus", NArm = TRUE)
s16sV1V3.4
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 551 taxa and 28 samples ]
    ## sample_data() Sample Data:       [ 28 samples by 4 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 551 taxa by 6 taxonomic ranks ]

``` r
plot_abundance = function(physeq,title = "",
                 Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("p__Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "sex",y = "Abundance",
                                 color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
                position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}

plotBefore = plot_abundance(s16sV1V3.4,"before")
```

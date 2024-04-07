read("ExperimentFunctions.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
sigstring = "0-2";


\\cpct ranges
b_ranges =
[
    1,12, 1;                    \\#best=3 previously up to 50
    1,12, 1;                    \\#best=3 previously up to 50
    1,12, 1;                    \\#best 3 previously up to 50
    5,20, 1;                   \\#best 13 up to 90: likely range [5, 25]
    7,22, 1;                   \\#best 18 max 90, likely range 5-30
    15,40, 2;                  \\#best 19 max 200. likely range [5, 25]
    50, 90, 3;\\10,25,2;        \\#best 75 max 250 likely range [90, 120]
    80, 120, 4;\\10,25,2;       \\#best 100 max 300 likely range [80, 100]
    20, 65, 3;\\10,25,2;        \\#best 40 max 300 likely range [15, 45]
    150, 250, 10; \\28, 48, 4;      \\#best 200 #max 1000, likely range [310, 370]
    60, 140, 10;  \\28, 48, 4;      \\#best 100 max 1000 likely range [200, 350]
    140, 240, 10; \\28, 48, 4;      \\#best 175 max 1000 likely range [50, 200]
    320, 450, 10;    \\#best 350 max 1000 likely range [1000, 1200]
    800, 1000, 25;    \\#best 900 2000 likely [1100]
    220, 400, 10;      \\#best 325 max1200 likely range [400, 1200]
    50,100,15;       \\#ceil(i/3)== 6
    50,100,15;       \\#ceil(i/3)== 6
    50,100,15       \\#ceil(i/3)== 6
];

/* log
b_ranges =
[
    5,50, 10;                    \\#previously up to 50
    1,50, 10;  \\1,7,1;          \\#previously up to 50
    1,50, 10;  \\1,7,1;          \\#previously up to 50
    2,25, 2;  \\3,10,1;         \\#up to 90: likely range [5, 25]
    2,25, 2;  \\3,10,1;         \\#max 90, likely range 5-30
    2,25, 2;  \\3,10,1;         \\#max 90. likely range [5, 25]
    90, 130, 5;\\10,25,2;       \\#max 250 likely range [90, 120]
    80, 100, 2;\\10,25,2;        \\#max 300 likely range [80, 100]
    15, 45, 3;\\10,25,2;        \\#max 300 likely range [15, 45]
    290, 370, 9; \\28, 48, 4;      \\#max 1500, likely range [310, 370]
    200, 350, 15; \\28, 48, 4;      \\#max1000 likely range [200, 350]
    50, 200, 50; \\28, 48, 4;      \\#max1000 likely range [50, 200]
    300, 600, 50;    \\#max4000 likely range [1000, 1200]
    1000, 1500, 50;    \\#2000 likely [1100]
    500, 700, 25;    \\#max1200 likely range [400, 1200]
    50,100,15;       \\#ceil(i/3)== 6
    50,100,15;       \\#ceil(i/3)== 6
    50,100,15       \\#ceil(i/3)== 6
];
*/
}

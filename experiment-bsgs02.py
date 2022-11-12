read("ExperimentFunctions.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
sigstring = "0-2";

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
    500, 1000, 50;    \\#max4000 likely range [1000, 1200]
    1000, 1500, 50;    \\#ceil(i/3)== 5
    500, 700, 25;    \\#max1200 likely range [400, 1200]
    50,100,15;       \\#ceil(i/3)== 6
    50,100,15;       \\#ceil(i/3)== 6
    50,100,15       \\#ceil(i/3)== 6
];
}


read("ExperimentFunctions.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
sigstring = "2-1";


b_ranges =
[
    16,25, 2;                        \\# [10, 50] best 18
    10, 22, 2;                       \\# [10, 50] best 17
    14, 24, 2;                       \\# [10, 50] best 20
    6, 16, 2;                        \\# [10, 50] best 11
    16, 26, 2;                       \\# [10, 50] best 22
    500,580, 10;                   \\# [100, 1000] best 520, terrible below 100
    230,300, 10;                   \\# [10, 1000] best 250
    250, 350, 10;                  \\# [10, 1000] best 275
    310, 380, 10;                  \\# [10, 1000] best 350
    110, 170,  10;                 \\# [10, 1000] best 130
    4300, 5000, 50;               \\#[1000, 6000] best 4600  sub 1400 is bad
    4300, 4900, 50;                 \\#[1000, 6000] best 4400
    4900, 5600, 50;                 \\#[1000, 6000] best 5200
    3100, 3500, 50;                 \\#[1000, 4000] best 3300
    4300, 4900, 50                  \\#[1000, 3000] best 4600
];

/*
b_ranges =
[
    1,12,1;          \\#100 -  [,3 ]
    1,12,1;          \\#100-  [,3 ]
    1,12,1;          \\#100-  [,3 ]
    5,25,2;          \\#70-  [15]
    5,25,2;         \\#110-  [13]
    5,30,2;          \\#90-  [11]
    20,60,3;         \\#400-  [30]
    30,85,5;      \\#1000-  [,60 ]
    20,100,8;         \\#700-  [, 40]
    190, 1000, 50;      \\#4000-  [,370 ]
    190, 1000, 50;      \\#4000-  [,200 ]
    50, 1000, 50;      \\#4000-  [, ]
    1000, 2500, 100;    \\#2000-  [600]
    1000, 2500, 100;    \\#2000-  [???, ]
    1000, 2000, 100;    \\#800 [300]
    100,600,50;       \\#ceil(i/3)== 6
    100,600,50;       \\#ceil(i/3)== 6
    100,600,50       \\#ceil(i/3)== 6
];
*/
}

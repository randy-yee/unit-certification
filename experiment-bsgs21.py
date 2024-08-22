
read("ExperimentFunctions.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
sigstring = "2-1";

b_ranges =
[
    16,25, 2;                        \\# [10, 50] best 18
    10, 22, 2;                       \\# [10, 50] best 18
    14, 24, 2;                       \\# [10, 50] best 20
    6, 16, 2;                        \\# [10, 50] best 10
    16, 26, 2;                       \\# [10, 50] best 22
    500,540, 5;                    \\# [100, 1000] best 520, terrible below 100
    240,320, 5;                   \\# [10, 1000] best 250
    280, 340, 5;                  \\# [10, 1000] best 300
    330, 370, 5;                  \\# [10, 1000] best 350
    100, 150,  5;                 \\# [10, 1000] best 120
    4300, 4800, 50;                 \\#[1000, 6000] best 4600  sub 1400 is bad
    4500, 4900, 50;                 \\#[1000, 6000] best 4800
    5000, 5800, 50;                 \\#[1000, 6000] best 5200
    3200, 3600, 50;                 \\#[1000, 4000] best 3450
    4300, 4800, 50;                  \\#[1000, 3000] best 4600
    10000, 20000, 2000            \\
];

b_ranges_round2 =
[
    70, 125, 5;                       \\# [10, 50] best 18
    120,170, 5;                       \\# [10, 50] best 18
    70, 110, 5;                       \\# [10, 50] best 20
    80, 120,5;                       \\# [10, 50] best 10
    50, 90, 5;                       \\# [10, 50] best 22
    650,900, 25;                   \\# [100, 1000] best 520, terrible below 100
    900,1150,25;                   \\# [10, 1000] best 250
    1100, 1500, 50;                \\# [10, 1000] best 300
    1700, 2000, 25;                \\# [10, 1000] best 350
    450, 650,  25;                 \\# [10, 1000] best 120
    30000, 37000, 1000;              \\#[1000, 6000] best 4600  sub 1400 is bad
    12000, 16000, 500;               \\#[1000, 6000] best 4800
    22000, 25000, 500;               \\#[1000, 6000] best 5200
    10000, 14000, 500;                  \\#[1000, 4000] best 3450
    20000, 25000, 500;                  \\#[1000, 3000] best 4600
    10000, 20000, 2000             \\
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

read("ExperimentFunctions.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
sigstring = "1-2";



b_ranges =
[
    3,7,  1;                      \\# [10, 40] best 4
    6,10, 1;                      \\# [10, 40] best 7
    4, 8, 1;                      \\# [10, 40] best 6
    5, 9, 1;                      \\# [10, 40] best 7
    3, 7, 1;                      \\# [10, 40] best 6
    40, 65, 2;                    \\# [10, 500] best 50
    75, 100, 2;                   \\# [10, 500] best 60
    45, 65,  2;                   \\# [10, 500] best 55
    80, 100, 2;                   \\# [10, 500] best 90
    75, 100, 2;                   \\# [10, 500] best 85
    1500, 1800, 25;               \\#[1000, 2000] best 1800
    1850, 2200, 25;                 \\#[1000, 2000] best 2050
    1100, 1400, 25;                 \\#[1000, 2000] best 1300
    1200, 1400, 25;                 \\#[1000, 2000] best 1400
    450, 750, 25;                 \\#[1000, 2000] best 600
    100,200,25;       \\#ceil(i/3)== 6
    100,200,25;       \\#ceil(i/3)== 6
    100,200,25       \\#ceil(i/3)== 6
];

/*
b_ranges =
[
    1,10,1;          \\#max70, likely [1..10]  4
    1,10,1;          \\#max90, likely [1..10]  3
    1,10,1;          \\#max90, likely [1..10]  3
    1,12,1;         \\#max100, likely [1..15]  4
    1,12,1;         \\#max100, likely [1..15]  7
    2,18,1;         \\#max100, likely [1..20]  7
    10,45,3;          \\#max150, likely [10..45] 25
    10,45,3;          \\#120, likely [1..20]
    3,30,2;           \\#max150, likely [10..60]
    90, 150, 6;        \\#max400, likely [10..60]
    40, 90, 5;          \\#max400, likely [60]
    60, 135, 5;         \\#max400, likely [105]
    600, 900, 30;    \\#ceil(i/3)== 5
    100, 300, 20;    \\#ceil(i/3)== 5
    100, 300, 20;    \\#ceil(i/3)== 5
    100,200,25;       \\#ceil(i/3)== 6
    100,200,25;       \\#ceil(i/3)== 6
    100,200,25       \\#ceil(i/3)== 6
];
*/
}

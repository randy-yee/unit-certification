read("ExperimentFunctions.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

{
sigstring = "3-0";



b_ranges =
[
    5,15,1;          \\#best 10 max 70, likely [1..10]  4
    4,15,1;          \\#best 10 max 90, likely [1..10]  3
    5,14,1;          \\#best 9  max 90, likely [1..10]  3
    25,35,1;         \\#best 31 max 200, likely [1..15]  4
    35,55,2;         \\#best 49 max 200, likely [1..15]  7
    22,38,2;         \\#best 32 max 200, likely [1..20]  7
    65,95,3;          \\ #best 74 max150, likely [10..45] 25
    65,105,4;          \\#best 93 max 180, likely [1..20]
    110,140,3;         \\#best 128 max 200-  [, 40]
    200, 260, 5;     \\##best 230  max 700 . likely range [200, 400]
    170, 290, 5;     \\#best 225 #max800 likely range [100, 300]
    130, 190, 5;      \\#best 160 #max700 likely range [1, 200]
    650, 850, 20;    \\#best 775 ceil(i/3)== 5
    1000, 1200, 20;    \\\best 1075, max 2000
    900, 1100, 20;    \\#best 950 max 1500# [1300]
    100,150,50;       \\#ceil(i/3)== 6
    100,150,50;       \\#ceil(i/3)== 6
    100,150,50       \\#ceil(i/3)== 6
];

/*
b_ranges =
[
    1,15,1;\\5,50,1;          \\#max 100. likely range [1, 15]
    1,11,1;\\5,50,1;          \\#max 100. likely range [1, 11]
    1,11,1;                   \\#max 100. likely range [1, 11]
    15,35,3;\\12,21,1;         \\#max 300. likely range [15,35]
    15,50,3;\\12,21,1;         \\\#max300 . likely range [15, 50]
    10,40,2;\\12,21,1;         \\#max300 . likely range [20]
    25,175,20;\\20,35,3;         \\#max300 . likely range [75 ]
    25,125,10;\\20,35,3;         \\#max . likely range [,75 ]
    30,170,10;\\20,35,3;         \\##max . likely range [,75 ]
    200, 380, 20;     \\##max700 . likely range [200, 400]
    100, 300, 20;     \\##max700 likely range [100, 300]
    60, 160, 10;      \\##max700 likely range [1, 200]
    600, 1100, 40;    \\#ceil(i/3)== 5
    1000, 1400, 50;    \\\#max2000# [1300]
    1000, 1400, 50;    \\#max2000# [1300]
    100,150,50;       \\#ceil(i/3)== 6
    100,150,50;       \\#ceil(i/3)== 6
    100,150,50       \\#ceil(i/3)== 6
];
*/
}

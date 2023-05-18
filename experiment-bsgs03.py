read("ExperimentFunctions.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Starting precision:
{
sigstring = "0-3";



b_ranges =
[
    1,15,1;          \\#best 4 max 90, likely [1..10]  4
    1,15,1;          \\#best 6 max90, likely [1..10]  3
    1,15,1;          \\#best 7 max90, likely [1..10]  3
    2,25,2;         \\#best 8 max 170, likely [1..15]  4
    2,20,2;         \\#best 6 max190, likely [1..15]  7
    5,35,3;         \\#best 10 max 200, likely [1..20]  7
    20,55,3;         \\#best 44 max150, likely [10..45] 25
    15,45,3;         \\#best 30 max 200, likely [1..20]
    5,45,3;         \\#best 23, max 400-  [, 40]
    115, 170, 5;     \\#best 140 max 300  likely range [,]
    50,  100, 5;     \\#best 60-75; max 1000  likely range [,]
    35,  95, 5;      \\#best 85; max 1000  likely range [,]
    410, 550, 10;     \\ #best 470 max1400 [900,1400]
    150, 250, 10;       \\#best 220 max 800 [300]
    250, 430, 10;      \\#best 340 max800 [300]
    25,100, 25;       \\#ceil(i/3)== 6
    25,100, 25;       \\#ceil(i/3)== 6
    25,100, 25       \\#ceil(i/3)== 6
];

/*
b_ranges =
[
    1,10, 1;           \\#max18  likely range [,]
    1,10, 1;           \\#max30  likely range [,]
    1,15, 1;         \\#max20  likely range [1,10]
    1,15, 2;         \\#max200  likely range [,]
    1,15, 2;         \\#max200  likely range [,]
    1,15, 2;         \\#max200  likely range [,]
    20,60, 4;         \\#max200  likely range [40,]
    10,40, 3;         \\#max117  likely range [20,]
    5,60, 3;         \\#max60  likely range [,10]
    100,  300, 20;     \\#max1000  likely range [,]
    10,  100, 10;     \\#max400  likely range [,]
    25,  200, 25;     \\#max400  likely range [,]
    900, 1500, 50;    \\#max1700 [1100]
    100, 600, 50;       \\#max800 [300]
    100, 600, 50;      \\#max800 [300]
    25,100, 25;       \\#ceil(i/3)== 6
    25,100, 25;       \\#ceil(i/3)== 6
    25,100, 25       \\#ceil(i/3)== 6
];
*/
}

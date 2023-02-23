read("ExperimentFunctions.py");
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Starting precision:
{
sigstring = "0-3";

b_ranges =
[
    1,18, 2;           \\#max18  likely range [,]
    1,25, 2;           \\#max30  likely range [,]
    1,20, 1;         \\#max20  likely range [1,10]
    5,60, 5;         \\#max200  likely range [,]
    5,60, 5;         \\#max200  likely range [,]
    5,60, 5;         \\#max200  likely range [,]
    20,60, 4;         \\#max200  likely range [40,]
    10,40, 3;         \\#max117  likely range [20,]
    5,60, 3;         \\#max60  likely range [,10]
    100,  300, 20;     \\#max1000  likely range [,]
    10,  100, 10;     \\#max400  likely range [,]
    25,  200, 25;     \\#max400  likely range [,]
    500, 1700, 200;       \\#ceil(i/3)== 5
    300, 1500, 100;       \\#ceil(i/3)== 5
    300, 1500, 100;      \\#ceil(i/3)== 5
    25,100, 25;       \\#ceil(i/3)== 6
    25,100, 25;       \\#ceil(i/3)== 6
    25,100, 25       \\#ceil(i/3)== 6
];
}

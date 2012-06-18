clear all; close all; clc

xp = [0 1 2 3];
x = 1.5;

L1 = ((x-xp(2))*(x-xp(3))*(x-xp(4)))/((xp(1)-xp(2))*(xp(1)-xp(3))*(xp(1)-xp(4)))
L2 = ((x-xp(1))*(x-xp(3))*(x-xp(4)))/((xp(2)-xp(1))*(xp(2)-xp(3))*(xp(2)-xp(4)))
L3 = ((x-xp(2))*(x-xp(1))*(x-xp(4)))/((xp(3)-xp(2))*(xp(3)-xp(1))*(xp(3)-xp(4)))
L4 = ((x-xp(2))*(x-xp(3))*(x-xp(1)))/((xp(4)-xp(2))*(xp(4)-xp(3))*(xp(4)-xp(1)))
function dummyOutput = savemsg(message)
    msg = msgbox(message, "Save notification");
    uiwait(msg);
    dummyOutput = NaN;
end
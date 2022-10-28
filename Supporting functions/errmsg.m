function dummyOutput = errmsg(message)
    msg = msgbox(message, "Error",'error');
    uiwait(msg);
    dummyOutput = NaN;
end
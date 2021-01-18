--This should be run with the command M2 --script m2r_server_script.m2 PORT#

if #commandLine < 4 then
  error("CommandLine was not given port.");

m2rintopenport = commandLine#3;

m2rintlog = (m2rinttoprint) -> (
  if m2rinttoprint == "" then (
    print("");
  ) else (
    m2rintdateinout = openInOut "!date \"+%Y-%m-%d %H:%M:%S: \"";
    m2rintthedate = read m2rintdateinout;
    print(concatenate(substring(m2rintthedate, 0, length m2rintthedate - 1), m2rinttoprint));
    close m2rintdateinout;
  );
);

m2rintlog("Session begins");

m2rintinout = openInOut concatenate("$:", m2rintopenport);
m2rintruncount = 0;
m2rintinout << "1.0.0" << "\n" << flush;

m2rintlog(concatenate("Connection received on port ", toString m2rintopenport));
m2rintlog("");

while true do (
  --while (not isReady m2rintinout) do (
  --  run "sleep 0.05";
  --);
  while (not isReady m2rintinout and not atEndOfFile m2rintinout) do (
    wait m2rintinout;
  );
  if (atEndOfFile m2rintinout) then break;
  m2rintinline = read m2rintinout;
  if (atEndOfFile m2rintinout or m2rintinline == "") then break;

  m2rintlog(concatenate("Command:\n", m2rintinline));

  m2rintretcode = 0;
  m2rintoutvalsucceeded = false;
  m2rintoutlinesucceeded = false;
  m2rintruncount = m2rintruncount + 1;

  try (
    m2rintoutval_m2rintruncount = value(m2rintinline);
    m2rintoutvalsucceeded = true;

    m2rintoutclass = class m2rintoutval_m2rintruncount;
    m2rintoutclassclass = class m2rintoutclass;

    m2rintvarname = "m2o" | toString(m2rintruncount);
    value(m2rintvarname | " = m2rintoutval_m2rintruncount;");

    m2rintoutline = toExternalString m2rintoutval_m2rintruncount;
    m2rintoutlinesucceeded = true;
  );

  if not m2rintoutvalsucceeded then (
    m2rintoutline = "Macaulay2 Error!";
    m2rintretcode = 1;
  ) else if not m2rintoutlinesucceeded then (
    m2rintoutline = "Macaulay2 toExternalString Error!";
    m2rintretcode = 2;
  );

  m2rintnumlines = 1 + #select("\n", m2rintoutline);
  m2rintoutinfo = concatenate(toString(m2rintretcode),
                          " ", toString(m2rintnumlines),
                          " ", toString(m2rintvarname),
                          " ", toString(m2rintoutclass),
                          " ", toString(m2rintoutclassclass));

  m2rintinout << m2rintoutinfo << "\n" << m2rintoutline << "\n" << flush;

  m2rintlog(concatenate("Command output:\n", m2rintoutinfo, "\n", m2rintoutline));
  m2rintlog("");
);

m2rintlog("Session ends");
close m2rintinout;

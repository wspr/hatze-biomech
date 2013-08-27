function person = person_generate(person,varargin)

ip = inputParser;
addParamValue(ip,'plot',false);
parse(ip,varargin{:})

person.plot = ip.Results.plot;

for ii = 1:person.N
  person = person.segment(ii).setup_fn(person,ii);
end

end
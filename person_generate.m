function person = person_generate(person,varargin)

ip = inputParser;
addParamValue(ip,'plot',false);
parse(ip,varargin{:})

person.plot = ip.Results.plot;

person = person_copy_angles(person);

for ii = 1:person.N
  person = person.segment(ii).setup_fn(person,ii);
end

end

function person = person_copy_angles(person)

person.segment( 1).angle = person.q(4:6);
person.segment( 2).angle = person.q(7:9);
person.segment( 3).angle = [person.q(10:11); 0];
person.segment( 4).angle = person.q(12:14);
person.segment( 5).angle = [person.q(15); 0; person.q(16)];
person.segment( 6).angle = [person.q(17:18); 0];
person.segment( 7).angle = [person.q(19:20); 0];
person.segment( 8).angle = person.q(21:23);
person.segment( 9).angle = [person.q(24); 0; person.q(25)];
person.segment(10).angle = [person.q(26:27); 0];
person.segment(11).angle = person.q(28:30);
person.segment(12).angle = person.q(31:33);
person.segment(13).angle = [person.q(34); 0; 0];
person.segment(14).angle = [person.q(35); 0; person.q(36)];
person.segment(15).angle = person.q(37:39);
person.segment(16).angle = [person.q(40); 0; 0];
person.segment(17).angle = [person.q(41); 0; person.q(42)];

for S = 1:person.N
  person.segment(S).Rlocal  = person.cardan_rotation(person.segment(S).angle);
end

end
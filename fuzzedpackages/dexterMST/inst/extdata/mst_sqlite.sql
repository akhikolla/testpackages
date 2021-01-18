pragma foreign_keys=1;--#split#--

CREATE TABLE Items(
	item_id VARCHAR(100) NOT NULL,
	PRIMARY KEY(item_id)
) 
WITHOUT ROWID;--#split#--

CREATE TABLE Scoring_rules(
	item_id VARCHAR(100) NOT NULL,
	response VARCHAR(100) NOT NULL,
	item_score INTEGER NOT NULL CHECK(typeof(item_score) = 'integer' AND item_score >= 0), 
	
	PRIMARY KEY (item_id, response),
	FOREIGN KEY (item_id) REFERENCES Items(item_id) ON UPDATE CASCADE ON DELETE CASCADE
) 
WITHOUT ROWID;--#split#--


CREATE TABLE Tests(
  test_id VARCHAR(100) NOT NULL PRIMARY KEY,
  routing VARCHAR(10) NOT NULL DEFAULT 'all' CHECK(routing IN('all','last') )
) 
WITHOUT ROWID;--#split#--

CREATE TABLE Booklets(
  test_id VARCHAR(100) NOT NULL,
	booklet_id VARCHAR(100) NOT NULL,
	
	PRIMARY KEY (test_id, booklet_id),
	FOREIGN KEY (test_id) REFERENCES Tests(test_id) ON UPDATE CASCADE ON DELETE CASCADE
) 
WITHOUT ROWID;--#split#-- 

CREATE TABLE Modules(
    test_id VARCHAR(100) NOT NULL,
	module_id VARCHAR(100) NOT NULL,
	
	PRIMARY KEY (test_id, module_id),
	FOREIGN KEY (test_id) REFERENCES Tests(test_id) ON UPDATE CASCADE ON DELETE CASCADE
) 
WITHOUT ROWID;--#split#--


CREATE TABLE Module_design(
	test_id VARCHAR(100) NOT NULL,
	module_id VARCHAR(100) NOT NULL,
	item_id VARCHAR(100) NOT NULL,
	item_position INTEGER NOT NULL CHECK(item_position >= 1),
	
	PRIMARY KEY (test_id, module_id, item_id),
	UNIQUE		(test_id, module_id, item_position),

	FOREIGN KEY (item_id) REFERENCES Items(item_id) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (test_id,  module_id) 
	  REFERENCES Modules(test_id,  module_id) ON UPDATE CASCADE ON DELETE CASCADE
) 
WITHOUT ROWID;--#split#--

CREATE TABLE Booklet_design(
  test_id VARCHAR(100) NOT NULL,
	booklet_id VARCHAR(100) NOT NULL,
	module_id VARCHAR(100) NOT NULL,
	module_nbr INTEGER NOT NULL CHECK(typeof(module_nbr) = 'integer' AND module_nbr >= 1),
	module_exit_score_min INTEGER NOT NULL 
	  CHECK(typeof(module_exit_score_min) = 'integer' AND module_exit_score_min >= 0),
	module_exit_score_max INTEGER NOT NULL 
	  CHECK(typeof(module_exit_score_max) = 'integer' AND module_exit_score_max >= module_exit_score_min),
	
	PRIMARY KEY (test_id, booklet_id, module_id ),
	UNIQUE      (test_id, booklet_id, module_nbr),
	
	FOREIGN KEY (test_id, module_id) REFERENCES Modules(test_id, module_id) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (test_id, booklet_id) REFERENCES Booklets(test_id, booklet_id) ON UPDATE CASCADE ON DELETE CASCADE
)
WITHOUT ROWID;--#split#--

-- needs many more elaborate contraints and triggers


 
CREATE TABLE Persons(
	person_id VARCHAR(100) NOT NULL,
	
	PRIMARY KEY(person_id)
);--#split#--



CREATE TABLE Administrations(
	person_id VARCHAR(100) NOT NULL,
	test_id VARCHAR(100) NOT NULL,
	booklet_id VARCHAR(100) NOT NULL,
	PRIMARY KEY (person_id, test_id),
	UNIQUE      (test_id, booklet_id, person_id)

	FOREIGN KEY (person_id) REFERENCES Persons(person_id) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (test_id, booklet_id) REFERENCES Booklets(test_id,booklet_id) ON UPDATE CASCADE ON DELETE CASCADE
) 
WITHOUT ROWID;--#split#--


CREATE TABLE Responses(
	person_id VARCHAR(100) NOT NULL,
	test_id VARCHAR(100) NOT NULL,
	booklet_id VARCHAR(100) NOT NULL,
	module_id VARCHAR(100) NOT NULL,
	item_id VARCHAR(100) NOT NULL,
	response VARCHAR(100) NOT NULL,
	
	PRIMARY KEY (test_id, person_id,  item_id),
	
	FOREIGN KEY (test_id, booklet_id, person_id) REFERENCES Administrations(test_id, booklet_id, person_id) ON UPDATE CASCADE ON DELETE CASCADE DEFERRABLE INITIALLY DEFERRED,
	FOREIGN KEY (test_id, booklet_id, module_id) REFERENCES Booklet_design(test_id, booklet_id, module_id) ON UPDATE CASCADE ON DELETE CASCADE DEFERRABLE INITIALLY DEFERRED,
	FOREIGN KEY (test_id, module_id, item_id) REFERENCES Module_design(test_id, module_id, item_id) ON UPDATE CASCADE ON DELETE CASCADE DEFERRABLE INITIALLY DEFERRED,
	FOREIGN KEY (item_id, response) REFERENCES Scoring_rules(item_id, response) ON UPDATE CASCADE ON DELETE CASCADE DEFERRABLE INITIALLY DEFERRED
	-- foreign key constraints deferred for speed 
) 
WITHOUT ROWID;

